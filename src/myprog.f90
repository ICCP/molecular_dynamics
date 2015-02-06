program ArgonGas

  use initialization
  use force
  use energy
  use const_temp

  implicit none

  integer, parameter :: particles = 864
  real(8), parameter :: density = 1.0
  real(8), parameter :: temp_target = 1.0

  integer, parameter :: boxes = nint((particles/4)**(1.0/3))
  real(8), parameter :: temperature = 1.0
  integer, parameter :: number_timesteps = 2000
  real(8), parameter :: time_step = 0.001
  real(8), parameter :: init_distance = (4.0/density)**(1.0/3)
  real(8), parameter :: length = boxes*init_distance
  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: ener_kin, ener_pot, temp_final
  integer :: i


 ! real(8) :: v_test(3)
    
  call initialize_position(position, boxes, particles, init_distance)
  call initialize_velocity(velocity, particles, temperature)
  call setting_cero_velocity(velocity, particles)
  call initialize_force(forces, particles)
  call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)

  open (unit=11,file='ener_kin_data.txt')
  open (unit=12,file='ener_pot_data.txt')
  open (unit=13,file='ener_tot_data.txt')
  !open (unit=14,file='temp_final_data.txt')

  do i = 0, number_timesteps
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     position(:,:) = position(:,:)+velocity(:,:)*time_step
     position(:,:) = modulo(position(:,:), length)
     call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     call calculate_kin_energy(velocity, particles, ener_kin)
     if (modulo(i,40) == 0 .and. i<700) then     
        call constant_temperature(particles, temp_target, velocity, temp_final)
     end if
     write(11,*) i, ener_kin
     write(12,*) i, ener_pot
     write(13,*) i, ener_kin+ener_pot
     !write(14,*) i, sum(velocity**2)/(3*(particles-1))
  end do
  
!  v_test(:)=0
!   do i=1, particles
!       v_test(:) = v_test(:)+velocity(:,i)
!    end do
!print*, v_test

end program ArgonGas
