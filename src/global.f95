module global

!Global Variable files

integer :: nxcells,nycells,nzcells    !Number of cells in x,y,z direction
integer :: ncells                     !Number of cells
integer :: ppc                        !Particles Per Cell
integer :: nprtl                      !Total Number of Particles
integer :: prtlnum                    !Particle Number
real(8) :: xcellscl,ycellscl,zcellscl !Width of fcc cell

!arrays that hold the positions and velocity of each particle
real(8), dimension(:,:,:),allocatable :: pos,vel,accel   
real(8), dimension(3) :: rand_vel = 0           !varible for assigning random velocities
real(8), dimension(3) :: prtl_accel = 0         !acceration on an individual particle

real(8) :: dt !timestep
integer :: NT !Number of timestep

!Array the hold the position of particles in the unit fcc cell
real(8), dimension(4,3) :: fcc


end module
