module initialization

  use global

  implicit none

  private open_files
  private read_parameters
  private allocate_initialize_var
  private initialize_position
  private initialize_velocity
  private setting_cero_velocity
  private gaussRandom
  private init_random_seed

  public Initialize_problem


contains

  subroutine Initialize_problem()


    call open_files()
    call read_parameters()
    call allocate_initialize_var()
    call initialize_position()
    call initialize_velocity()
    call setting_cero_velocity()
    pair_corr = 0._8
    pressures = 0._8
    pot_energy = 0._8
    vel_sqr = 0._8
    flag = 1

!!$    call initplot ('lightblue', 800,800, 'out.ps', 1)
!!$    call Framing (0._8, 0._8, length, length)
!!$    call putstopbutton()
  
  end subroutine Initialize_problem

  subroutine open_files()
    
    open (unit=11,file='ener_kin_data.txt', status='REPLACE')
    open (unit=12,file='ener_pot_data.txt', status='REPLACE')
    open (unit=13,file='ener_tot_data.txt', status='REPLACE')
    open (unit=14,file='temp_final_data.txt', status='REPLACE')
    open (unit=15,file='pair_corre_data.txt', status='REPLACE')
    open (unit=16,file='pressure_data.txt', status='REPLACE')
    open (unit=17,file='Parameters.txt', status='OLD')
    
  end subroutine open_files

  subroutine read_parameters()

    read(17,*) num_particles
    read(17,*) density
    read(17,*) temp_target
    read(17,*) number_timesteps

    print*, num_particles, density, temp_target, number_timesteps

  end subroutine read_parameters

  subroutine allocate_initialize_var()

    if(.not.allocated(position))allocate(position(3,num_particles))
    if(.not.allocated(velocity))allocate(velocity(3,num_particles))
    if(.not.allocated(forces))allocate(forces(3,num_particles))

    boxes = nint((num_particles/4)**(1.0/3))
    length = boxes*(4.0/density)**(1.0/3)

    if(.not.allocated(pair_corr))allocate(pair_corr(nint(length/2*100)))
    if(.not.allocated(pressures))allocate(pressures(number_timesteps))
    if(.not.allocated(pot_energy))allocate(pot_energy(number_timesteps))
    if(.not.allocated(vel_sqr))allocate(vel_sqr(number_timesteps))


  end subroutine allocate_initialize_var
 
  subroutine initialize_position()
    
    real(8) :: init_distance
    integer :: i, j, k, n=1

    init_distance = length/boxes

    do i = 0, boxes-1
       do j = 0, boxes-1
          do k = 0, boxes-1
             !First particle in origin
             position(1,n) = i*init_distance
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance
             n = n+1
             !Second particle in face k=1
             position(1,n) = i*init_distance+init_distance/2.0
             position(2,n) = j*init_distance+init_distance/2.0
             position(3,n) = k*init_distance
             n = n+1
             !Third in face i=1
             position(1,n) = i*init_distance            
             position(2,n) = j*init_distance+init_distance/2.0
             position(3,n) = k*init_distance+init_distance/2.0
             n = n+1
             !Fourth in face j=1
             position(1,n) = i*init_distance+init_distance/2.0    
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance+init_distance/2.0
             n = n+1
          end do
       end do
    end do
   
  end subroutine initialize_position
   
  subroutine initialize_velocity()

    integer :: i,j

    call init_random_seed


    do j = 1, num_particles
       do i = 1, 3
          velocity(i,j) = gaussRandom()
       end do
    end do

  end subroutine initialize_velocity

  real(8) function gaussRandom() result (rand_vel)
    
    real(8), parameter :: Temperature = 1.0
    real(8) :: random1, random2, probability

    call random_number(random1)
    call random_number(random2)

    ! Calculate a random number between -3T/2 and 3T/2 and its probability
    random1 = random1*3.0*sqrt(Temperature) - 3.0*sqrt(Temperature)/2.0
    probability = (1/sqrt(2.0*Temperature*PI))*exp(-random1**2/(2.0*Temperature))

    !if the second number is lower than the probability 
    !you accept it, if not you do it again
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

  subroutine setting_cero_velocity()

    real(8) :: v_tot(3)
    integer :: j

    v_tot = 0

    do j=1, num_particles
       v_tot = v_tot + velocity(:,j)
    end do

    do j=1, num_particles
       velocity(:,j)= velocity(:,j) - v_tot/num_particles
    end do

  end subroutine setting_cero_velocity

end module initialization
