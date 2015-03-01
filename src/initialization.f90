module initialization

  use global

  implicit none

  private initialize_position
  private initialize_velocity
  private setting_cero_velocity
  private gaussRandom
  private init_random_seed

  public Initialize_problem
  public open_files
  public read_parameters
  public allocate_initialize_var

contains

  subroutine Initialize_problem()

    !$$ To initialize the problem we do the following:

    !$$ Fisrt we open all needed files where we are going to
    !$$ store all the physical quantities
    call open_files()
    
    !$$ Second we need to read the parameters the user can change
    !$$ and initialize some variables
    call read_parameters()
    call allocate_initialize_var()
    pair_corr = 0._8
    pressures = 0._8
    pot_energy = 0._8
    vel_sqr = 0._8
    flag = 1

    !$$ Third for the system we initialize the positions,
    !$$ and the random the random velocities 
    call initialize_position()
    call initialize_velocity()
    call setting_cero_velocity()

!!$    call initplot ('lightblue', 800,800, 'out.ps', 1)
!!$    call Framing (0._8, 0._8, length, length)
!!$    call putstopbutton()
  
  end subroutine Initialize_problem

  subroutine open_files()
    
    open (unit=15,file='pair_corre_data.txt', status='REPLACE')
    open (unit=17,file='Parameters.txt', status='OLD')
    
  end subroutine open_files

  subroutine read_parameters()

    read(17,*) num_particles
    read(17,*) density
    read(17,*) temp_target
    read(17,*) number_timesteps
    read(17,*) r_cutoff

  end subroutine read_parameters

  subroutine allocate_initialize_var()

    !$$ Here we allocate the dimension of our arrays
    !$$ in dependance with the parameters the user gave
    if(.not.allocated(position))allocate(position(3,num_particles))
    if(.not.allocated(velocity))allocate(velocity(3,num_particles))
    if(.not.allocated(forces))allocate(forces(3,num_particles))

    boxes = nint((num_particles/4)**(1.0/3))
    length = boxes*(4.0/density)**(1.0/3)
    delta_r = length/400

    if(.not.allocated(pair_corr))allocate(pair_corr(nint(length/(2*delta_r))))
    if(.not.allocated(pressures))allocate(pressures(number_timesteps))
    if(.not.allocated(pot_energy))allocate(pot_energy(number_timesteps))
    if(.not.allocated(vel_sqr))allocate(vel_sqr(number_timesteps))


  end subroutine allocate_initialize_var
 
  subroutine initialize_position()
    
    !$$ In this subroutine we initialize the position of
    !$$ the particles. We set them in a FCC latice.

    real(8) :: init_distance
    integer :: i, j, k, n=1

    !$$ This is the length of one side of the FCC
    init_distance = length/boxes

    do i = 0, boxes-1
       do j = 0, boxes-1
          do k = 0, boxes-1
             !$$ First particle in origin
             position(1,n) = i*init_distance
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance
             n = n+1
             !$$ Second particle in face k=1
             position(1,n) = i*init_distance+init_distance/2.0
             position(2,n) = j*init_distance+init_distance/2.0
             position(3,n) = k*init_distance
             n = n+1
             !$$ Third in face i=1
             position(1,n) = i*init_distance            
             position(2,n) = j*init_distance+init_distance/2.0
             position(3,n) = k*init_distance+init_distance/2.0
             n = n+1
             !$$ Fourth in face j=1
             position(1,n) = i*init_distance+init_distance/2.0    
             position(2,n) = j*init_distance
             position(3,n) = k*init_distance+init_distance/2.0
             n = n+1
             !$$ We repeat this procidure through all the boxes
          end do
       end do
    end do
   
  end subroutine initialize_position
   
  subroutine initialize_velocity()

    !$$ In this subroutine we assign a random number to 
    !$$ every veolcity in 3D of each particle

    integer :: i,j

    call init_random_seed

    do j = 1, num_particles
       do i = 1, 3
          velocity(i,j) = gaussRandom()
       end do
    end do

  end subroutine initialize_velocity

  real(8) function gaussRandom() result (rand_vel)
    
    !$$ In this function we calculate a random number
    !$$ following a Gaussian Distribution

    real(8), parameter :: Temperature = 1.0
    real(8) :: random1, random2, probability

    !$$ Create 2 random numbers and safe them in random1 and random2
    call random_number(random1)
    call random_number(random2)

    !$$ Calculate a random number between -3T/2 and 3T/2 and its probability
    random1 = random1*3.0*sqrt(Temperature) - 3.0*sqrt(Temperature)/2.0
    probability = (1/sqrt(2.0*Temperature*PI))*exp(-random1**2/(2.0*Temperature))

    !$$ If the second number is lower than the probability of the first one 
    !$$ you accept it, if not you get other numbers and do it again
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

    v_tot = 0._8

    !$$ First we calculate the total velocity of the system
    do j=1, num_particles
       v_tot = v_tot + velocity(:,j)
    end do

    !$$ Then we just substract this velocity to get 0 total velocity
    do j=1, num_particles
       velocity(:,j)= velocity(:,j) - v_tot/num_particles
    end do

  end subroutine setting_cero_velocity

end module initialization
