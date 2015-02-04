program ArgonGas

  use initialization
  use force

  implicit none

  integer, parameter :: boxes = 2, particles = 32
  integer, parameter :: temperature = 3
  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)

  call initialize_position(position, boxes)
  call initialize_velocity(velocity, particles, temperature)
  call setting_cero_velocity(velocity, particles)
  call calculate_force(forces, particles, position, boxes)
      
  write(*,"(3E10.2)") position !3 reals of width 10 and 2 decimals
  print *, "fuerzas" 
  write(*,"(3E10.2)") forces


end program
