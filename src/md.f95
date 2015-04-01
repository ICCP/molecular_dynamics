program md

  implicit none
  
  integer :: nprtl          !Number of particle in a direction
  integer :: ncells         !Number of cells in one dimension
  integer :: NT
  real(8) :: boxlength      !Length of one box
  real(8) :: neighbor_dist  !Distance to nearest neighbor
  real(8) :: bounds
  real(8),dimension(:,:),allocatable :: pos,vel,accel
  
  NT = 1000      !Number of Time Steps
  ncells = 1
  nprtl  = 4*ncells**3
  neighbor_dist = 2.d0**(1.d0/6.d0)-0.1
  boxlength = 2.d0 / sqrt(2.d0) * neighbor_dist
  bounds = ncells*boxlength

  allocate(pos(nprtl,3))
  allocate(vel(nprtl,3))
  allocate(accel(nprtl,3))
  
  pos = 0
  vel = 0 
  accel = 0

  open(unit = 2, file = 'pos.xyz', status = 'unknown')
 
  print*,'neighbor_dist',neighbor_dist
  print*,'boxlength', boxlength 
  print*,'bounds:',bounds
  
  call bld_lattice(pos,ncells,boxlength,nprtl)
  call xyz_file_printout(pos,nprtl,1)
  !call force_calc(pos,nprtl,accel)
  !call xyz_file_printout(accel,nprtl,2)
  call verlet_velocity_integration(pos,vel,accel,nprtl,NT,bounds)

end program 

subroutine bld_lattice(pos_prtl,numcells_wide,cell_width,numprtl) !{{{
  
  implicit none

  !Input
  real(8) :: cell_width
  integer :: numcells_wide
  integer :: numprtl
  real(8),dimension(numprtl,3) :: pos_prtl

  !Internal
  integer :: ii,jj,kk,ll,prtl_count
  real(8),dimension(4,3) :: fcc


  !FCC Cordinates
  fcc(1,:) = [0.d0,0.d0,0.d0]
  fcc(2,:) = [0.d0,cell_width/2.d0,cell_width/2.d0]
  fcc(3,:) = [cell_width/2.d0,0.d0,cell_width/2.d0]
  fcc(4,:) = [cell_width/2.d0,cell_width/2.d0,0.d0]
 
  prtl_count = 0

  do ii = 1,numcells_wide
      do jj = 1,numcells_wide
          do kk = 1, numcells_wide
              do ll = 1,4
                  prtl_count = (ii-1)*numcells_wide**3 + (jj-1)*numcells_wide**2 + (kk-1)*numcells_wide + ll
                  pos_prtl(prtl_count,1) = fcc(ll,1) + (ii-1)*cell_width
                  pos_prtl(prtl_count,2) = fcc(ll,2) + (jj-1)*cell_width
                  pos_prtl(prtl_count,3) = fcc(ll,3) + (kk-1)*cell_width
              end do 
          end do
      end do 
  end do 


end subroutine !}}}

subroutine xyz_file_printout(pos_prtl,numprtl,timestep) !{{{
  
  implicit none

  !Input Variables
  integer :: numprtl
  integer :: timestep
  real(8),dimension(numprtl,3) :: pos_prtl
  
  !Internal Variables
  integer :: ii


  write(2,*),numprtl
  write(2,*),'Particles at time step:',timestep
  !Internal Variables
  do ii = 1,numprtl
      write(2,*),ii,pos_prtl(ii,:)
  end do 

end subroutine !}}}

subroutine force_calc(pos_prtl,numprtl,boundry,force)!{{{

  implicit none

  !Input Variables
  integer :: numprtl
  real(8),dimension(numprtl,3) :: pos_prtl,force
  real(8) :: boundry
  !Internal Variables
  integer :: ii, jj 
  real(8),dimension(3) :: dr
  real(8) :: drsq       !value of 
  real(8) :: force_mag
  real(8) :: min_dist
  min_dist = 10000
  do ii = 1,numprtl-1
      do jj = ii+1,numprtl
          dr = pos_prtl(jj,:) - pos_prtl(ii,:)
          dr = dr - NINT(dr/boundry)*boundry
          drsq = dot_product(dr,dr)
          !print out the smallest interations distance
          if (min_dist > sqrt(drsq)) min_dist = sqrt(drsq)
          
          force_mag = (-48.d0/drsq**7 + 24/drsq**4)
          force(ii,:) = force(ii,:) + force_mag*dr
          force(jj,:) = force(jj,:) - force_mag*dr
      end do 
  end do 

print*,'min_dist:', min_dist
end subroutine!}}}

subroutine verlet_velocity_integration(pos_prtl,vel_prtl,accel_prtl,numprtl,numts,boundry)!{{{

  implicit none

  !Input Variables
  integer :: numprtl
  integer :: numts
  real(8),dimension(numprtl,3) :: pos_prtl,vel_prtl,accel_prtl
  real(8) :: boundry
  !Internal Variables
  integer :: ii
  real(8) :: dt = 0.004   !Not a paramerter I want to change
  real(8),dimension(numprtl,3) :: vel_p,vel_pp,accel_p,accel_pp,pos_p,pos_pp
  
  vel_p    = 0.d0
  vel_pp   = 0.d0
  accel_p  = 0.d0
  accel_pp = 0.d0
  pos_p    = 0.d0
  pos_pp   = 0.d0


  do ii = 1,numts

      !Calculate vel for t+.5
      !vel_prtl = vel_prtl + 0.5*accel_prtl*dt                
      
      !Calculate pos for t+1
      !pos_prtl = pos_prtl + vel_prtl*dt                      
      
      !Boundry Condition on position
      !pos_prtl = modulo(pos_prtl,boundry)
      !Past time steps move to past past time steps
      pos_pp   = pos_p
      vel_pp   = vel_p 
      accel_pp = accel_p
      
      !Current time steps go to past time steps
      pos_p    = pos_prtl
      vel_p    = vel_prtl
      accel_p  = accel_prtl

      !Calculate vel for t+.5
      vel_prtl = vel_p + 0.5*accel_p*dt                
      
      !Calculate pos for t+1
      pos_prtl = pos_p + vel_p*dt 
      pos_prtl = pos_prtl + 0.5d0*dt*dt*accel_p

      accel_prtl = 0.d0
      
      !Recalculate Force
      call force_calc(pos_prtl,numprtl,boundry,accel_prtl)   !Calculate acceleration  for t+1
      
      !Update Velocity
      vel_prtl = vel_prtl + 0.5*accel_prtl*dt                !Calculate vel for t+1

      call xyz_file_printout(pos_prtl,numprtl,ii)
  end do 
end subroutine!}}}

subroutine verlet_velocity(pos_prtl,vel_prtl,accel_prtl,numprtl,numts,boundry)!{{{

  implicit none

  !Input Variables
  integer :: numprtl
  integer :: numts
  real(8),dimension(numprtl,3) :: pos_prtl,vel_prtl,accel_prtl
  real(8) :: boundry
  
  !Internal Variables
  integer :: ii,jj,kk
  real(8) :: dt = 0.004   !Not a paramerter I want to change
  real(8),dimension(3) :: dr
  real(8) :: drsq
  real(8) :: force_mag
  real(8) :: min_dist
  real(8),dimension(numprtl,3) :: vel_p,vel_pp,accel_p,accel_pp,pos_p,pos_pp
  
  vel_p    = 0.d0
  vel_pp   = 0.d0
  accel_p  = 0.d0
  accel_pp = 0.d0
  pos_p    = 0.d0
  pos_pp   = 0.d0
  
  do ii = 1,numts
      !Past time steps move to past past time steps
      pos_pp   = pos_p
      vel_pp   = vel_p 
      accel_pp = accel_p
      
      !Current time steps go to past time steps
      pos_p    = pos_prtl
      vel_p    = vel_prtl
      accel_p  = accel_prtl

      !Calculate vel for t+.5
      vel_prtl = vel_p + 0.5*accel_p*dt                
      
      !Calculate pos for t+1
      pos_prtl = pos_p + vel_p*dt 
      pos_prtl = pos_prtl + 0.5d0*dt*dt*accel_p
      
      !Boundry Condition on position
      pos_prtl = modulo(pos_prtl,boundry)
      
      !Recalculate Force
      accel_prtl = 0.d0
      min_dist = 10000
      do jj = 1,numprtl-1
          do kk = jj+1,numprtl
              dr = pos_prtl(kk,:) - pos_prtl(jj,:)
              dr = dr - dble(NINT(dr/boundry)*boundry)
              drsq = dot_product(dr,dr)
              !print out the smallest interations distance
              if (min_dist > sqrt(drsq)) min_dist = sqrt(drsq)
          
              force_mag = (-48.d0/drsq**(7.d0) + 24.d0/drsq**(4.d0))
              accel_prtl(jj,:) = accel_prtl(jj,:) + force_mag*dr
              accel_prtl(kk,:) = accel_prtl(kk,:) - force_mag*dr
          end do 
      end do 

      print*,'min_dist:', min_dist
      !Update Velocity
      vel_prtl = vel_prtl + 0.5*accel_prtl*dt                !Calculate vel for t+1

      call xyz_file_printout(pos_prtl,numprtl,ii)
  end do 


end subroutine!}}}

