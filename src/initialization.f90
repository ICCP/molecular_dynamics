module initialization

  implicit none

  private gaussRandom
  private init_random_seed

  public initialize_position
  public initialize_velocity
  public setting_cero_velocity

contains
  
  subroutine initialize_position(position, L, num_part, init_distance)
    
    integer, intent(in) :: L, num_part
    real(8), intent(in) :: init_distance
    real(8), intent(inout) :: position(3,num_part)
    integer :: i, j, k, n=1
    
    do i = 0, L-1
       do j = 0, L-1
          do k = 0, L-1
             position(1,n) = i*init_distance            !First particle in origin
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance
             n = n+1
             position(1,n) = i*init_distance+2.0**(-1.0/3)!Second particle in face k=1
             position(2,n) = j*init_distance+2.0**(-1.0/3)
             position(3,n) = k*init_distance
             n = n+1
             position(1,n) = i*init_distance            !Third in face i=1
             position(2,n) = j*init_distance+2.0**(-1.0/3)
             position(3,n) = k*init_distance+2.0**(-1.0/3)
             n = n+1
             position(1,n) = i*init_distance+2.0**(-1.0/3)    !Fourth in face j=1
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance+2.0**(-1.0/3)
             n = n+1
          end do
       end do
    end do
   
  end subroutine initialize_position
   
  subroutine initialize_velocity(velocity, N, T)

    real(8), intent(in) :: T
    integer, intent(in) :: N
    real(8), intent(inout) :: velocity(3,N)

    integer :: i,j

    call init_random_seed

    do i = 1, 3
       do j = 1, N
          velocity(i,j) = gaussRandom(T)
       end do
    end do

  end subroutine initialize_velocity

  real(8) function gaussRandom(Temperature) result(rand_vel)
    
    real(8), intent(in) :: Temperature
    real(8), parameter :: PI = 4*atan(1.0)
    real(8) :: random1, random2, probability

    call random_number(random1)
    call random_number(random2)

    random1 = random1*3.0*sqrt(Temperature) - 3.0*sqrt(Temperature)/2.0
    probability = (1/sqrt(2.0*Temperature*PI))*exp(-random1**2/(2.0*Temperature))

    do while (probability < random2)
       call random_number(random1)
       call random_number(random2)
       
       random1 = random1*3.0*sqrt(Temperature) - 3.0*sqrt(Temperature)/2.0
       probability = (1/sqrt(2.0*Temperature*PI))*exp(-random1**2/(2.0*Temperature))
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

  subroutine setting_cero_velocity(velocity, N)

    integer , intent(in) :: N
    real(8), intent(inout) :: velocity(3,N)
    real(8) :: v_tot(3)
    integer :: j

    v_tot(:)=0

    do j=1, N
       v_tot(:) = v_tot(:)+velocity(:,j)
    end do

    do j=1, N
       velocity(:,j)= velocity(:,j)-v_tot(:)/N
    end do

  end subroutine setting_cero_velocity

end module initialization
