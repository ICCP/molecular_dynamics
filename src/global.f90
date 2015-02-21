module global

  implicit none

  integer :: num_particles = 864
  real(8) :: density = 0.88
  real(8) :: temp_target = 1

  integer :: number_timesteps = 2000
  real(8), parameter :: time_step = 0.004
  real(8), parameter :: r_cutoff = 3.2
  integer, parameter :: time_cut = 500

  real(8), allocatable :: position(:,:)
  real(8), allocatable :: velocity(:,:)
  real(8), allocatable :: forces(:,:)

  real(8), parameter :: PI = 4*atan(1.0)
  integer :: boxes
  real(8) :: length
  integer :: flag

  real(8), allocatable :: pair_corr(:)
  real(8), allocatable :: pressures(:)

  real(8), allocatable :: pot_energy(:)
  real(8), allocatable :: vel_sqr(:)
  
contains 


end module global
