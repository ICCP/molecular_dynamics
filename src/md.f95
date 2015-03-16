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
  fcc(3,:) = (/0.5,0.5,0.0/)
  fcc(4,:) = (/0.5,0.0,0.5/)
 
  NT = 1

  allocate(pos(nprtl,3,1))   !allocate position
  allocate(vel(nprtl,3,1))   !allocate velocity
  allocate(accel(nprtl,3,1)) !allocate acceration

  call bld_lattice  !Create inital position of gas particles
  call writepos     !write positions to files
  call writevel     !write velocities to files



  call force_lj(fcc(1,:),fcc(2,:),prtl_force_lj)

  print*,'prtl_force_lj',prtl_force_lj
  print*,'Norm force', norm2(prtl_force_lj)
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
        pos(prtlnum,1,1) = fcc(ll,1) + (ii-1)*xcellscl
        pos(prtlnum,2,1) = fcc(ll,2) + (jj-1)*ycellscl
        pos(prtlnum,3,1) = fcc(ll,3) + (kk-1)*zcellscl
  
        !Create and scale random velocities
        call random_number(rand_vel)
        call scalerand(rand_vel)

        !Assign random velocities
        vel(prtlnum,1:3,1) = rand_vel
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
    write (1,*),ii,pos(ii,:,1)
  end do
end subroutine !}}}

subroutine writevel!{{{
  !Functionality - Write Velocity out to file
  
  use global
   
  implicit none

  open (unit = 2, file = 'vel.out', status = 'unknown')
  do ii = 1,nprtl
    write (2,*),ii,vel(ii,:,1)
  end do

end subroutine !}}}

subroutine scalerand(randvel)  !{{{

  use global 

  implicit none

  real(8),dimension(3) :: randvel

  randvel = 2*randvel-1
  
end subroutine  !}}}

subroutine force_lj(pos1,pos2,force) !{{{
  
  !Calculates the force due to the Lenard Jones Potiential 
  !Input: pos1,pos2- position of the first and second particle
  
  use global

  implicit none
  
  !input variabels 
  real(8),dimension(3),intent(in) :: pos1,pos2  !position of particle 1 and 2
  real(8),dimension(3),intent(out) :: force      !force between particles
  
  !Internal Variable
  real(8),dimension(3) :: r ,r_norm  !Distance Between Radius
  real(8) :: force_mag               !Magnitude of the Force

  r = pos1 - pos2             !Finding the Distance Between the Points
  r_norm = r/norm2(r)         !Normalized Distance

  !Lenard Jones Potential 
  force_mag = 24.d0*((2.d0/(norm2(r))**13)+(-1.d0/(norm2(r))**7))

  !Direction Force
  force = r_norm * force_mag

end subroutine !}}}

subroutine verlet_integration(pos,vel,accel,NT)
  !------------------------------------------------------------------------
  !Function - Calculate the postion of all the particles after NT timesteps
  !
  !Input:  pos - position of every particle over NT timesteps
  !        vel - velocity of every particle over NT timesteps
  !        accel - acceration of every particle over NT timesteps
  !        NT - number of timesteps 
  !Output: pos
  !------------------------------------------------------------------------
  
  real(8), dimension(:,:,:) :: pos, vel, accel
  integer, intent(in) :: NT 

  pos(:,:,2) = pos(:,:,1) + vel(:,:,1) * dt + 0.5*accel(:,:,1) * dt**2

  do ii = 2 , NT-1
    pos(:,:,ii+1) = 2.d0*pos(:,:,ii) - pos(:,:,ii-1) + accel(:,:,ii)*dt**2
  end do 

end subroutine











