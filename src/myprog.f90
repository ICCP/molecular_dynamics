program ArgonGas

  use initialization
  use force
  use energy
  use const_temp

  implicit none

  integer, parameter :: boxes = 6, particles = 864
  real(8), parameter :: temperature = 0.1
  real(8), parameter :: total_time = 3
  real(8), parameter :: time_step = 0.001
  real(8), parameter :: init_distance = 2.0**(2.0/3)
  real(8), parameter :: length = boxes*init_distance
  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: ener_kin, ener_pot, temp_target = 1.0
  integer :: totaltimeint
  integer :: i
 ! real(8) :: v_test(3)
  totaltimeint = nint(total_time/time_step)

  call initialize_position(position, boxes, particles, init_distance)
  call initialize_velocity(velocity, particles, temperature)
  call setting_cero_velocity(velocity, particles)
  call initialize_force(forces, particles)
     call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)

  open (unit=7,file='ener_kin_data.txt')
  open (unit=8,file='ener_pot_data.txt')
  open (unit=9,file='ener_tot_data.txt')

  do i = 0, totaltimeint
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     position(:,:) = position(:,:)+velocity(:,:)*time_step
     position(:,:) = modulo(position(:,:), length)
     call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     call calculate_kin_energy(velocity, particles, ener_kin)
     if (i>500 .and. i<1500) then     
        call constant_temperature(particles, temp_target, velocity)
     end if
     write(7,*) i, ener_kin
     write(8,*) i, ener_pot
     write(9,*) i, ener_kin+ener_pot
  end do

!  v_test(:)=0
!   do i=1, particles
!       v_test(:) = v_test(:)+velocity(:,i)
!    end do
!print*, v_test

end program ArgonGas
