program ArgonGas

  use initialization
  use time_evol
  use force
  use energy
  use pressure
  use pair_correl

  implicit none

  integer, parameter :: num_particles = 2048
  real(8), parameter :: density = 1.2
  real(8), parameter :: temp_target = 0.5

  integer, parameter :: number_timesteps = 2000

  real(8) :: position(3,num_particles)
  real(8) :: velocity(3,num_particles)
  real(8) :: forces(3,num_particles)

  call Initialize_problem(num_particles, density, position, velocity, forces) 
    
  call Time_evolution(num_particles, density, temp_target, number_timesteps, position, velocity, forces)!, pair_corre, pressure)

!!$  call EndPlot()

end program ArgonGas
