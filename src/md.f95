program md

  implicit none
  
  integer :: nxcells,nycells,nzcells,ncells,ppc,nprtl,prtlnum
  real(8),dimension(:,:),allocatable :: pos,vel
  real(8),dimension(4,3) :: fcc
  integer :: ii,jj,kk
  !Program Settings
  nxcells = 10 !Number of cells in the X directions
  nycells = 10 !Number of cells in the Y directions
  nzcells = 10 !Number of cells in the Z direction
  ncells = nxcells*nycells*nzcells !Number of total boxes
  ppc = 4      !particle per cell
  nprtl = ncells*ppc !number of particles
  fcc(1,:) = (/0.0,0.0,0.0/)
  fcc(2,:) = (/0.0,0.5,0.5/)
  fcc(3,:) = (/0.5,0.5,0.0 /)
  fcc(4,:) = (/0.5,0.0,0.5 /)
  allocate(pos(nprtl,3))

  do ii = 1,ncells
    do jj = 1,ppc
        pos(ii,1) = fcc(jj,1) + mod(ii,nxcells)
        pos(ii,2) = fcc(jj,2) + 
        pos(ii,3) = fcc(jj,3) + 
    end do
  end do 

end program
