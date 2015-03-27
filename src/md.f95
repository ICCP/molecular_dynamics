program md

  use global

  implicit none
  
  !Program Settings
  nxcells = 1    !Number of cells in the X direction
  nycells = 1    !Number of cells in the Y direction
  nzcells = 1    !Number of cells in the Z direction
  
  xcellscl = 1.0 !Width of cell in X direction
  ycellscl = 1.0 !Width of cell in y direction
  zcellscl = 1.0 !Width of cell in z direction
  
  scalefactor = 1.12  

  ncells = nxcells*nycells*nzcells !Number of Total Boxes
  ppc = 4                          !Particle per Cell
  
  !Set up boundry conditions
  xbound = nxcells*xcellscl*scalefactor 
  ybound = nycells*ycellscl*scalefactor
  zbound = nzcells*zcellscl*scalefactor
  
  !nprtl = ncells*ppc !number of particles
  nprtl = 1

  dt = 0.004*1.d0
  
  !FCC Cordinates
  fcc(1,:) = (/0.0,0.0,0.0/)
  fcc(2,:) = (/0.0,0.5,0.5/)
  fcc(3,:) = (/0.5,0.5,0.0/)
  fcc(4,:) = (/0.5,0.0,0.5/)
 
  NT = 10000  !Number of timesteps

  allocate(pos(nprtl,3,NT))   !allocate position
  allocate(vel(nprtl,3,NT))   !allocate velocity
  allocate(accel(nprtl,3,NT)) !allocate acceration

  
  !-----------Main Program-----------!
  call bld_lattice_two_prtl
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

  integer :: ii, jj, kk, ll

  do ii = 1,nzcells
    do jj = 1,nycells
      do kk = 1,nxcells
        do ll = 1,ppc
        prtlnum = (ii-1)*nycells*nxcells*ppc+(jj-1)*nxcells*ppc+(kk-1)*ppc + ll
        
        !Set inital position
        pos(prtlnum,1,1) = (fcc(ll,1) + (ii-1)*xcellscl)*scalefactor
        pos(prtlnum,2,1) = (fcc(ll,2) + (jj-1)*ycellscl)*scalefactor
        pos(prtlnum,3,1) = (fcc(ll,3) + (kk-1)*zcellscl)*scalefactor
  
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
  !pos(2,:,1) = [0.d0,0.d0,1.2*1.d0]

  vel(1,:,1) = [0.d0,0.d0,.5*1.d0]
  !vel(2,:,1) = [0.d0,0.d0,0.d0]
end subroutine !}}}

subroutine writepos !{{{
  !Functionality - write position out to file

  use global
  
  implicit none

  integer :: ii,jj

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

  integer :: ii,jj

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

  randvel = randvel
  
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
  r(1) = r(1) - NINT(r(1)/xbound)*xbound
  r(2) = r(2) - NINT(r(2)/ybound)*ybound
  r(3) = r(3) - NINT(r(3)/zbound)*zbound
  
  !Lenard Jones force 
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

  implicit none

  integer :: ii

  call accel_calc(1)
  !first time step iteration

  pos(:,:,2) = pos(:,:,1) + vel(:,:,1) * dt + 0.5*accel(:,:,1)*dt**2
  pos(:,1,2) = modulo(pos(:,1,2),xbound)
  pos(:,2,2) = modulo(pos(:,2,2),ybound)
  pos(:,3,2) = modulo(pos(:,3,2),zbound)
  print*,'pos1',pos(1,:,1) 

  do ii = 2 , NT
      call accel_calc(ii)
      pos(:,:,ii+1) = 2.d0*pos(:,:,ii) - pos(:,:,ii-1) + accel(:,:,ii)*dt**2
      pos(:,1,ii+1) = modulo(pos(:,1,ii+1),xbound)
      pos(:,2,ii+1) = modulo(pos(:,2,ii+1),ybound)
      pos(:,3,ii+1) = modulo(pos(:,3,ii+1),zbound)
      print*,"pos1",pos(1,:,ii)
  end do 

end subroutine !}}}

subroutine accel_calc(it) !{{{
  !Function: Calculates the acceleration of every particle from the iteration of every other particle
  !
  !input: it - current time step indexer
  !
  !Global Variables - 

  use global
 
  implicit none

  !internal variable
  real(8), dimension(3) :: prtl_force_lj !particle force from Lennard-Jones
  integer :: it                          !current iteration
  integer :: ii, jj


  do ii = 1, nprtl
      do jj = 1, nprtl 
          if (ii .ne. jj) then  
              call force_lj(pos(ii,:,it),pos(jj,:,it),prtl_force_lj) 
              accel(ii,:,it) = accel(ii,:,it) + prtl_force_lj
          end if 
      end do 
  end do 


end subroutine !}}}

subroutine sys_energy(it) !{{{
  !Function: Calculates the system energy 
  !
  !input: it - current time step ind
  !
  !Global Variables - 

  use global
 
  implicit none

  !internal variable
  real(8), Lenard_Jones_Potiential                              !Potiential Energy 
  integer :: it                          !current iteration
  integer :: ii, jj


  do ii = 1, nprtl
      do jj = 1, nprtl 
          if (ii .ne. jj) then  
              call force_lj(pos(ii,:,it),pos(jj,:,it),prtl_force_lj) 
              accel(ii,:,it) = accel(ii,:,it) + prtl_force_lj
          end if 
      end do 
  end do 



    0.d5*

end subroutine!}}}

subroutine U_lj(pos1,pos2,p_energy)!{{{

  !Finding the distance between points and adding periodic boundry conditions
  r = pos1 - pos2
  r(1) = r(1) - NINT(r(1)/xbound)*xbound
  r(2) = r(2) - NINT(r(2)/ybound)*ybound
  r(3) = r(3) - NINT(r(3)/zbound)*zbound

  !Lenard Jones Potiential 
  p_energy = 4.d0 * ((dot_product(r,r))**-6 - (dot_product(r,r))**-3)

end subroutine U_lj!}}}




