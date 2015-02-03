program ArgonGas

  implicit none

  integer, parameter :: boxes = 2, particles = 32
  integer, parameter :: temperature = 3
  real(8) :: position(3,32), velocity(3,32)


  call initialize_position(position, boxes)
  call initialize_velocity(velocity, particles, temperature)
      
  write(*,"(3E10.2)") position !3 reals of width 10 and 2 decimals
  print *, "Velocidades" 
  write(*,"(3E10.2)") velocity

contains

  subroutine initialize_position(position, L)

    integer, intent(in) :: L
    real(8), intent(inout) :: position(3,32)
    integer :: i, j, k, n=1

    do i = 0, L-1
       do j = 0, L-1
          do k = 0, L-1
             position(1,n) = i*2**(2.0/3)            !First particle in origin
             position(2,n) = j*2**(2.0/3)
             position(3,n) = k*2**(2.0/3)
             n = n+1
             position(1,n) = i*2**(2.0/3)+2**(-1.0/3)!Second particle in face k=1
             position(2,n) = j*2**(2.0/3)+2**(-1.0/3)
             position(3,n) = k*2**(2.0/3)
             n = n+1
             position(1,n) = i*2**(2.0/3)            !Third in face i=1
             position(2,n) = j*2**(2.0/3)+2**(-1.0/3)
             position(3,n) = k*2**(2.0/3)+2**(-1.0/3)
             n = n+1
             position(1,n) = i*2**(2.0/3)+2**(-1.0/3)    !Fourth in face j=1
             position(2,n) = j*2**(2.0/3)
             position(3,n) = k*2**(2.0/3)+2**(-1.0/3)
             n = n+1
          end do
       end do
    end do
   
  end subroutine initialize_position
   
  subroutine initialize_velocity(velocity, N, T)

    real(8), intent(inout) :: velocity(3,32)
    integer, intent(in) :: N, T

    integer :: i,j,k =8888

    do i = 1, 3
       do j = 1, N
          k = k+100
          velocity(i,j) = gaussRandom(T,k )
       end do
    end do

  end subroutine initialize_velocity

  real(8) function gaussRandom(Temperature, seed) result(rand_vel)
    
    integer, intent(in) :: Temperature, seed
    real(8), parameter :: PI=4*atan(1.0)
    real(8) :: random1, random2, probability

    call init_random_seed
    call random_number(random1)
    call random_number(random2)

    random1 = random1*6*sqrt(3.0*Temperature) - 3*sqrt(3.0*Temperature)
    probability = (1/sqrt(6.0*Temperature*PI))*exp(-random1**2/6.0*Temperature)

    print *, random1, probability, random2

    do while (probability < random2)
       call random_number(random2)
    end do

    rand_vel = random1


  end function gaussRandom

  subroutine init_random_seed()
    implicit none

    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)

    allocate(seed(n))
    open(newunit=un, file="/dev/urandom", access="stream",&
         form="unformatted", action="read", status="old", &
         iostat=istat)

    if (istat == 0) then
       read(un) seed
       close(un)

    else
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)

          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = ieor(s, pid)

       if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
    end if

    call random_seed(put=seed)

  end subroutine init_random_seed

end program
