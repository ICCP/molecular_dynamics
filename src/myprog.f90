program ArgonGas

  use initialization
  use force
  use energy
  use const_temp

  implicit none

  integer, parameter :: particles = 864
  real(8), parameter :: density = 1
  real(8), parameter :: temp_target = 1

  integer, parameter :: boxes = nint((particles/4)**(1.0/3))
  integer, parameter :: number_timesteps = 2000
  real(8), parameter :: time_step = 0.004
  real(8), parameter :: length = boxes*(4.0/density)**(1.0/3)

  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: pair_corre(nint(length*sqrt(3.0)*0.5/0.01))
  real(8) :: ener_kin, ener_pot, virial_function, pressure

  integer :: i,j, flag

  flag = 1
    
  call initialize_position(position, particles, density)
  call initialize_velocity(velocity, particles)
  call setting_cero_velocity(velocity, particles)
  forces = 0._8

!!$  call initplot ('lightblue', 800,800, 'out.ps', 1)
!!$  call Framing (0._8, 0._8, length, length)
!!$  call putstopbutton()

  call calculate_force(forces, particles, position, length, ener_pot, pair_corre, flag, virial_function)

  open (unit=11,file='ener_kin_data.txt')
  open (unit=12,file='ener_pot_data.txt')
  open (unit=13,file='ener_tot_data.txt')
  open (unit=14,file='temp_final_data.txt')
  open (unit=15,file='pair_corre_data.txt')
  open (unit=16,file='pressure_data.txt')

  do i = 0, number_timesteps
     velocity = velocity + forces*time_step/2.0
     position = modulo(position+velocity*time_step, length)
     call calculate_force(forces, particles, position, length, ener_pot, pair_corre, flag, virial_function)
     velocity = velocity + forces*time_step/2.0
!!$     call plot_particles(position, particles)
     ener_kin = sum(velocity**2*0.5)
     call thermostat(particles, temp_target, velocity, i)
     pressure = particles*temp_target/(length**3) + (virial_function/(3*particles*(length**3)))

     if (i == 700) then
        flag = 0
     end if

     write(11,*) i, ener_kin
     write(12,*) i, ener_pot
     write(13,*) i, ener_kin+ener_pot
     write(16,*) i, pressure
  end do
 
  do j = 1, nint(length*sqrt(3.0)*0.5/0.01)
     write(15,*) j/100.0, pair_corre(j)*2.0*(length**3)/(particles*(particles-1)*(number_timesteps-700))
  end do
  
!!$  !Proove to velocities equal to cero
!!$  v_test(:)=0
!!$  do i=1, particles
!!$     v_test(:) = v_test(:)+velocity(:,i)
!!$  end do
!!$  print*, v_test

!!$  call EndPlot()

!!$contains
!!$  
!!$  subroutine plot_particles(pos, n)
!!$    integer :: i, n
!!$    real(8) :: pos(3,n)
!!$    
!!$    do i= 1, n
!!$       call setpoint(pos(1,i), pos(2,i))
!!$    end do
!!$  end subroutine plot_particles

end program ArgonGas
