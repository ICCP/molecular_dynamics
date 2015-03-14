program md

  use global

  implicit none
  
  !Program Settings
  nxcells = 10   !Number of cells in the X directions
  nycells = 10   !Number of cells in the Y directions
  nzcells = 10   !Number of cells in the Z direction
  xcellscl = 1.0 !Width of cell in X direction
  ycellscl = 1.0 !Width of cell in y direction
  zcellscl = 1.0 !Width of cell in z direction
  ncells = nxcells*nycells*nzcells !Number of total boxes
  ppc = 4            !particle per cell
  nprtl = ncells*ppc !number of particles

  !FCC Cordinates
  fcc(1,:) = (/0.0,0.0,0.0/)
  fcc(2,:) = (/0.0,0.5,0.5/)
  fcc(3,:) = (/0.5,0.5,0.0 /)
  fcc(4,:) = (/0.5,0.0,0.5 /)
  
  allocate(pos(nprtl,3))   !allocate position
  allocate(vel(nprtl,3))   !allocate velocity
  allocate(accel(nprtl,3)) !allocate acceration

  call bld_lattice  !Create inital position of gas particles
  call writepos     !write positions to files
  call writevel     !write velocities to files


end program

subroutine bld_lattice !{{{
  !Functionality: Creates the inital position of all the particles in the box.
  !
  !Input: pos - position array that hold the cartesian cordinates of the particles.
  !             Array on input should be zeros.
  !Output: pos -
  
  use global

  implicit none

  do ii = 1,nzcells
    do jj = 1,nycells
      do kk = 1,nxcells
        do ll = 1,ppc
        prtlnum = (ii-1)*nycells*nxcells*ppc+(jj-1)*nxcells*ppc+(kk-1)*ppc + ll
        
        !Set inital position
        pos(prtlnum,1) = fcc(ll,1) + (ii-1)*xcellscl
        pos(prtlnum,2) = fcc(ll,2) + (jj-1)*ycellscl
        pos(prtlnum,3) = fcc(ll,3) + (kk-1)*zcellscl
  
        !Create and scale random velocities
        call random_number(rand_vel)
        call scalerand(rand_vel)

        !Assign random velocities
        vel(prtlnum,1:3) = rand_vel
        end do 
      end do 
    end do
  end do

end subroutine !}}}

subroutine writepos !{{{
  !Functionality - write position out to file

  use global
  
  implicit none

  open (unit = 1, file = 'pos.out', status = 'unknown')
  do ii = 1,nprtl
    write (1,*),ii,pos(ii,:)
  end do
end subroutine !}}}

subroutine writevel!{{{
  !Functionality - Write Velocity out to file
  
  use global
   
  implicit none

  open (unit = 2, file = 'vel.out', status = 'unknown')
  do ii = 1,nprtl
    write (2,*),ii,vel(ii,:)
  end do

end subroutine !}}}

subroutine scalerand(randvel)  !{{{

  use global 

  implicit none

  real(8),dimension(3) :: randvel

  randvel = 2*randvel-1
  
end subroutine  !}}}
