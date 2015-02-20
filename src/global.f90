module global

  implicit none

  integer, parameter :: num_particles = 864
  real(8), parameter :: density = 0.88
  real(8), parameter :: temp_target = 1.095

  integer, parameter :: number_timesteps = 2000

  real(8) :: position(3,num_particles)
  real(8) :: velocity(3,num_particles)
  real(8) :: forces(3,num_particles)

  real(8), parameter :: PI = 4*atan(1.0)

  real(8), parameter :: time_step = 0.004

  integer :: boxes
  real(8) :: length
  
contains 


end module global
