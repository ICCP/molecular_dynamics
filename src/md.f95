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
  !nprtl = ncells*ppc !number of particles
  nprtl = 2
  dt = 1.d0
  !FCC Cordinates
  fcc(1,:) = (/0.0,0.0,0.0/)
  fcc(2,:) = (/0.0,0.5,0.5/)
  fcc(3,:) = (/0.5,0.5,0.0/)
  fcc(4,:) = (/0.5,0.0,0.5/)
 
  NT = 3

  allocate(pos(nprtl,3,NT))   !allocate position
  allocate(vel(nprtl,3,NT))   !allocate velocity
  allocate(accel(nprtl,3,NT)) !allocate acceration

  call bld_lattice_two_prtl
  print*,'verlet_integration'
  call verlet_integration
  call writepos
  !call bld_lattice  !Create inital position of gas particles
  !call force_lj(fcc(2,:),fcc(1,:),force_test)
  !print*,'force_test:',force_test
  !call accel_calc(1)



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
        pos(prtlnum,1,1) = (fcc(ll,1) + (ii-1)*xcellscl)*1.90637
        pos(prtlnum,2,1) = (fcc(ll,2) + (jj-1)*ycellscl)*1.90637
        pos(prtlnum,3,1) = (fcc(ll,3) + (kk-1)*zcellscl)*1.90637
  
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

subroutine bld_lattice_two_prtl !{{{

  use global

  implicit none

  pos(1,:,1) = [0.d0,0.d0,0.d0]    
  pos(2,:,1) = [0.d0,0.d0,1.90637*1.d0]

  vel(1,:,1) = [0.d0,0.d0,-1.d0]
  vel(2,:,1) = [0.d0,0.d0,1.d0]
end subroutine !}}}

subroutine writepos !{{{
  !Functionality - write position out to file

  use global
  
  implicit none

  open (unit = 1, file = 'pos.out', status = 'unknown')
  do ii = 1,nprtl
      do jj = 1, NT
          write (1,*),ii,pos(ii,:,jj)
      end do
  end do
end subroutine !}}}

subroutine writevel!{{{
  !Functionality - Write Velocity out to file
  
  use global
   
  implicit none

  open (unit = 2, file = 'vel.out', status = 'unknown')
  do ii = 1,nprtl
      do jj = 1,NT
          write (2,*),ii,vel(ii,:,jj)
      end do
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
  real(8),dimension(3) :: r  !Distance Between Radius
  real(8) :: force_mag       !Magnitude of the Force

  r = pos1 - pos2             !Finding the Distance Between the Points

  !Lenard Jones Potential 
  force_mag = 24.d0*((2.d0/(dot_product(r,r))**7)+(-1.d0/(dot_product(r,r))**4))

  !Direction Force
  force = r * force_mag

end subroutine !}}}

subroutine verlet_integration !{{{
  !------------------------------------------------------------------------
  !Function - Calculate the postion of all the particles after NT timesteps
  !
  !Input:  pos - position of every particle over NT timesteps
  !        vel - velocity of every particle over NT timesteps
  !        accel - acceration of every particle over NT timesteps
  !        NT - number of timesteps 
  !Output: pos
  !------------------------------------------------------------------------
  
  use global 
  print*,'First iteration'
  call accel_calc(1)
  print*,'accel time = 1 :',accel(:,:,1)
  pos(:,:,2) = pos(:,:,1) + vel(:,:,1) * dt + 0.5*accel(:,:,1) * dt**2
  do ii = 2 , NT
      print*,'iteration',ii
      call accel_calc(ii)
      pos(:,:,ii+1) = 2.d0*pos(:,:,ii) - pos(:,:,ii-1) + accel(:,:,ii)*dt**2
  end do 

end subroutine !}}}

subroutine accel_calc(it) !{{{

  use global
  
  !internal variable
  real(8), dimension(3) :: prtl_force_lj !particle force from lenard jones
  integer :: it                          !current iteration
 
 do ii = 1, nprtl
      do jj = 1, nprtl 
          if (ii .ne. jj) then  
              print*,'ii,jj', ii,jj
              call force_lj(pos(ii,:,it),pos(jj,:,it),prtl_force_lj) 
              print*,'prtl_force', prtl_force_lj
              accel(ii,:,it) = accel(ii,:,it) + prtl_force_lj
          end if 
      end do 
  end do 


end subroutine !}}}









