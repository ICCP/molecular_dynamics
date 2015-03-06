program md

  use global

  implicit none
  

  !integer :: nxcells,nycells,nzcells,ncells,ppc,nprtl,prtlnum
  !real(8),dimension(:,:),allocatable :: pos,vel
  !real(8),dimension(4,3) :: fcc
  !integer :: ii,jj,kk,ll
  !real(8) :: xcellscl,ycellscl,zcellscl

  !Program Settings
  nxcells = 10   !Number of cells in the X directions
  nycells = 10   !Number of cells in the Y directions
  nzcells = 10   !Number of cells in the Z direction
  xcellscl = 1.0 !Width of cell in X direction
  ycellscl = 1.0 !Width of cell in y direction
  zcellscl = 1.0 !Width of cell in z direction
  ncells = nxcells*nycells*nzcells !Number of total boxes
  ppc = 4      !particle per cell
  nprtl = ncells*ppc !number of particles

  !FCC Cordinates
  fcc(1,:) = (/0.0,0.0,0.0/)
  fcc(2,:) = (/0.0,0.5,0.5/)
  fcc(3,:) = (/0.5,0.5,0.0 /)
  fcc(4,:) = (/0.5,0.0,0.5 /)
  allocate(pos(nprtl,3))
  
  !create lattice
  do ii = 1,nzcells
    do jj = 1,nycells
      do kk = 1,nxcells
        do ll = 1,ppc
        prtlnum = (ii-1)*nycells*nxcells*ppc+(jj-1)*nxcells*ppc+(kk-1)*ppc + ll
        pos(prtlnum,1) = fcc(ll,1) + (ii-1)*xcellscl
        pos(prtlnum,2) = fcc(ll,2) + (jj-1)*ycellscl
        pos(prtlnum,3) = fcc(ll,3) + (kk-1)*zcellscl
        end do 
      end do 
    end do
  end do 

  !write lattice to file
  open (unit = 1, file = 'pos.out', status = 'unknown')
  do ii = 1,nprtl
    write (1,*),ii,pos(ii,:)
  end do


end program

subroutine bld_lattice
  !Functionality: Creates the inital position of all the particles in the box.
  !
  !Input: pos - position array that hold the cartesian cordinates of the particles.
  !             Array on input should be zeros.
  !Output: pos -
  
  implicit none

end subroutine
