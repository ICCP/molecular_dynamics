program ArgonGas

  use initialization
  use force
  use energy

  implicit none

  integer, parameter :: boxes = 2, particles = 32
  real(8), parameter :: temperature = 1d0
  real(8), parameter :: total_time = 1
  real(8), parameter :: time_step = 0.001
  real(8), parameter :: length = boxes*(2.0**(2.0/3))
  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: ener_kin, ener_pot

  integer :: totaltimeint
  integer :: i!, j, k

  totaltimeint = nint(total_time/time_step)

  call initialize_position(position, boxes, particles)
  call initialize_velocity(velocity, particles, temperature)
  !velocity = 0d0
  call setting_cero_velocity(velocity, particles)

  call initialize_force(forces, particles)
  call calculate_force(forces, particles, position, boxes, ener_pot)
  call calculate_kin_energy(velocity, particles, ener_kin)
  print*, ener_kin, ener_pot, ener_kin+ener_pot

  !print *, sqrt(sum((position(:,1)-position(:,5))**2))

  !print *, "velocidades" 
  !write(*,"(3ES10.2)") forces
  !write(*,*) length

  do i = 0, totaltimeint
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     position(:,:) = position(:,:)+velocity(:,:)*time_step
     position(:,:) = modulo(position(:,:), length)
     call calculate_force(forces, particles, position, boxes, ener_pot)
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     call calculate_kin_energy(velocity, particles, ener_kin)
     print*, i, ener_kin, ener_pot, ener_kin+ener_pot
  end do

      
  !print *, "posiciones"
  !write(*,"(3E10.2)") position !3 reals of width 10 and 2 decimals
  !print *, "velocidades" 
  !write(*,"(3E10.2)") velocity
  !print *, "fuerzas" 
  !write(*,"(3E10.2)") forces


end program
