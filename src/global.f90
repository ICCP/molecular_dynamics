module global

  implicit none

  !$$ Global variables to use through all the program

  !$$ This are the variables provided by the user
  integer :: num_particles
  real(8) :: density
  real(8) :: temp_target
  integer :: number_timesteps

  !$$ This are the variables defining the time evolution
  !$$ and cut off of the potential 
  real(8), parameter :: time_step = 0.004
  real(8) :: r_cutoff
  integer, parameter :: time_cut = 500

  !$$ This are the basic arrays for the evolution of the problem
  real(8), allocatable :: position(:,:)
  real(8), allocatable :: velocity(:,:)
  real(8), allocatable :: forces(:,:)

  !$$ Just some variables that we are going 
  !$$ to use through all the program
  real(8), parameter :: PI = 4*atan(1.0)
  integer :: boxes
  real(8) :: length
  real(8) :: delta_r
  integer :: flag

  !$$ Physical quantities in an array
  real(8), allocatable :: pair_corr(:)
  real(8), allocatable :: pressures(:)
  real(8), allocatable :: pot_energy(:)
  real(8), allocatable :: vel_sqr(:)

end module global
