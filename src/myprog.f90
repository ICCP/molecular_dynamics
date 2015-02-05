program ArgonGas

  use initialization
  use force
  use energy
  use const_temp

  implicit none

  integer, parameter :: boxes = 2, particles = 32
  real(8), parameter :: temperature = 1d0
  real(8), parameter :: temp_target = 1.5
  real(8), parameter :: total_time = 1
  real(8), parameter :: time_step = 0.001
  real(8), parameter :: init_distance = 2.**(2./3)
  real(8), parameter :: length = boxes*init_distance
  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: ener_kin, ener_pot, temp_final = 0.0

  integer :: totaltimeint
  integer :: i!, j, k

  totaltimeint = nint(total_time/time_step)

  call initialize_position(position, boxes, particles, init_distance)
  call initialize_velocity(velocity, particles, temperature)
  !velocity = 0d0
  call setting_cero_velocity(velocity, particles)

  call initialize_force(forces, particles)
  call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)
  call calculate_kin_energy(velocity, particles, ener_kin)

  open (unit=7,file='ener_kin_data.txt')
  open (unit=8,file='ener_pot_data.txt')
  open (unit=9,file='ener_tot_data.txt')
  open (unit=10,file='temp_data.txt')

  do i = 0, totaltimeint
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     position(:,:) = position(:,:)+velocity(:,:)*time_step
     position(:,:) = modulo(position(:,:), length)
     call calculate_force(forces, particles, position, boxes, ener_pot, init_distance)
     velocity(:,:) = velocity(:,:) + forces(:,:)*time_step/2.0
     call calculate_kin_energy(velocity, particles, ener_kin)
     call constant_temperature(particles, temp_target, velocity, temp_final)
     write(7,*) i, ener_kin
     write(8,*) i, ener_pot
     write(9,*) i, ener_kin+ener_pot
     write(10,*) i, temp_final
  end do

  print*, ener_kin, ener_pot, ener_kin+ener_pot, temp_final

end program
