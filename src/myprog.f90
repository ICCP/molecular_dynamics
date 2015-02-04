program ArgonGas

 ! use md_plot
  use initialization
  use force

  implicit none

  integer, parameter :: boxes = 2, particles = 32, mass = 1
  integer, parameter :: temperature = 3
  real(8) :: positions(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: boxSize = boxes*2**(2.0/3), timestep = 1

  call initialize_position(positions, boxes)
  call initialize_velocity(velocity, particles, temperature)
  call setting_cero_velocity(velocity, particles)
  call calculate_force(forces, particles, positions, boxes)
      
  write(*,"(3E10.2)") positions !3 reals of width 10 and 2 decimals
  print *, "fuerzas" 
  write(*,"(3E10.2)") forces

  call md_plot_init(boxSize)
  call md_plot_points(positions)
  
  call md_plot_close()

end program
