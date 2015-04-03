program md

  implicit none
  
  integer :: nprtl          !Number of particle in a direction
  integer :: ncells         !Number of cells in one dimension
  integer :: NT
  real(8) :: boxlength      !Length of one box
  real(8) :: neighbor_dist  !Distance to nearest neighbor
  real(8) :: bounds
  real(8) :: temperature
  real(8),dimension(:,:),allocatable :: pos,vel,accel
  print*,'here'

  temperature = 2.d0
  NT = 500     !Number of Time Steps
  ncells = 8
  nprtl  = 4*ncells**3
  neighbor_dist = 2.d0**(1.d0/6.d0)
  boxlength = 2.d0 / sqrt(2.d0) * neighbor_dist
  bounds = ncells*boxlength

  allocate(pos(nprtl,3))
  allocate(vel(nprtl,3))
  allocate(accel(nprtl,3))
  
  pos = 0
  vel = 0 
  accel = 0

  open(unit = 50, file = 'pos.xyz', status = 'unknown')
  open(unit = 51, file = 'energy.out', status = 'unknown')
  open(unit = 52, file = 'pressure.out', status = 'unknown')
  open(unit = 53, file = 'temperature.out',status = 'unknown')
  open(unit = 54, file = 'raddist.out', status = 'unknown')

  print*,'neighbor_dist',neighbor_dist
  print*,'boxlength', boxlength 
  print*,'bounds:',bounds
  print*,'bould lattice'
  call bld_lattice(pos,vel,ncells,boxlength,nprtl,temperature)
  call xyz_file_printout(pos,nprtl,1)
  call verlet_velocity_integration(pos,vel,accel,nprtl,NT,bounds,temperature)

end program 

subroutine bld_lattice(pos_prtl,vel_prtl,numcells_wide,cell_width,numprtl,temp) !{{{
  
  implicit none

  !Input
  real(8) :: cell_width
  real(8) :: temp
  integer :: numcells_wide
  integer :: numprtl
  real(8),dimension(numprtl,3) :: pos_prtl,vel_prtl

  !Internal
  integer :: ii,jj,kk,ll,prtl_count
  real(8),dimension(4,3) :: fcc
  real(8) :: urv1,urv2,urv3,urv4 !four uniform random values
  real(8) :: nrv1,nrv2,nrv3,nrv4 !four normal random values
  real(8),dimension(3) :: vel_sum

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
                  prtl_count = prtl_count + 1
                  pos_prtl(prtl_count,1) = fcc(ll,1) + (ii-1)*cell_width
                  pos_prtl(prtl_count,2) = fcc(ll,2) + (jj-1)*cell_width
                  pos_prtl(prtl_count,3) = fcc(ll,3) + (kk-1)*cell_width
                  
                  !generate random values for velocity
                  call random_number(urv1)
                  call random_number(urv2)
                  call random_number(urv3)
                  call random_number(urv4)

                  nrv1 = sqrt(-2.d0*log(urv1))*cos(6.28318530718d0*urv2)
                  nrv2 = sqrt(-2.d0*log(urv2))*sin(6.28318530718d0*urv1)
                  nrv3 = sqrt(-2.d0*log(urv3))*cos(6.28318530718d0*urv4)
                  nrv4 = sqrt(-2.d0*log(urv4))*sin(6.28318530718d0*urv3)

                  vel_prtl(prtl_count,1) = nrv1
                  vel_prtl(prtl_count,2) = nrv2
                  vel_prtl(prtl_count,3) = nrv3

              end do 
          end do
      end do 
  end do 

  vel_sum(1) = sum(vel_prtl(:,1))
  vel_sum(2) = sum(vel_prtl(:,2))
  vel_sum(3) = sum(vel_prtl(:,3))

  vel_prtl(:,1) = vel_prtl(:,1) - vel_sum(1)/numprtl
  vel_prtl(:,2) = vel_prtl(:,2) - vel_sum(2)/numprtl
  vel_prtl(:,3) = vel_prtl(:,3) - vel_sum(3)/numprtl
  call scale_vel(vel_prtl,numprtl,0.5d0*sum(vel_prtl**2),temp)

end subroutine !}}}

subroutine xyz_file_printout(pos_prtl,numprtl,timestep) !{{{
  
  implicit none

  !Input Variables
  integer :: numprtl
  integer :: timestep
  real(8),dimension(numprtl,3) :: pos_prtl
  
  !Internal Variables
  integer :: ii


  write(50,*),numprtl
  write(50,*),'Particles at time step:',timestep
  !Internal Variables
  do ii = 1,numprtl
      write(50,*),ii,pos_prtl(ii,:)
  end do 

end subroutine !}}}

subroutine force_calc(pos_prtl,numprtl,boundry,force,ljp,press_temp)!{{{

  implicit none

  !Input Variables
  integer :: numprtl
  real(8),dimension(numprtl,3) :: pos_prtl,force
  real(8) :: boundry
  real(8) :: ljp
  real(8),dimension(100) :: rad_dist
  real(8) :: press_temp

  !Internal Variables
  integer :: ii, jj ,kk
  real(8),dimension(3) :: dr
  real(8) :: drsq       !value of 
  real(8) :: drsqrt
  real(8) :: force_mag
  real(8) :: min_dist
  press_temp = 0.d0
  min_dist = 10000
  ljp = 0.d0
  rad_dist = 0.d0
  do ii = 1,numprtl-1
      do jj = ii+1,numprtl
          dr = pos_prtl(jj,:) - pos_prtl(ii,:)
          dr = dr - NINT(dr/boundry)*boundry
          drsq = dot_product(dr,dr)
          drsqrt = sqrt(drsq)
          !print*,'drsqrt:',drsqrt
          !print out the smallest interations distance
          if (min_dist > drsqrt) min_dist = drsqrt
          
          !Calculating Lennard Jones Potiential
          ljp = ljp + 4.d0 * ( drsq**(-6.d0) - drsq**(-3.d0))
          
          !foce calculation
          force_mag = (-48.d0/drsq**7 + 24/drsq**4)
          force(ii,:) = force(ii,:) + force_mag*dr
          force(jj,:) = force(jj,:) - force_mag*dr
          
          !Pressure calculation
          press_temp = press_temp + force_mag*dot_product(dr,dr)

          !Radial Distibution
          do kk = 1,100
              !print*,'Lower Bounds:',((kk-1)*(boundry/10.d0)+0.7d0) 
              !print*,'Upper Bounds:',((kk)*(boundry/10.d0)+0.7d0)
              if((((kk-1)*(boundry/100.d0)+0.7d0) .lt. drsqrt) .and. (((kk)*(boundry/100.d0)+0.7d0) .gt. drsqrt)) then 
                
                  rad_dist(kk) = rad_dist(kk) + 1
              end if 
          end do 

      end do 
  end do 

  do kk = 1,100
      write(54,*),kk,rad_dist(kk)
  end do 

print*,'min_dist:', min_dist
end subroutine!}}}

subroutine verlet_velocity_integration(pos_prtl,vel_prtl,accel_prtl,numprtl,numts,boundry,temp)!{{{

  implicit none

  !Input Variables
  integer :: numprtl
  integer :: numts
  real(8),dimension(numprtl,3) :: pos_prtl,vel_prtl,accel_prtl
  real(8) :: boundry
  real(8) :: temp
  
  !Internal Variables
  integer :: ii
  real(8) :: dt = 0.004   !Not a paramerter I want to change
  real(8),dimension(numprtl,3) :: vel_p,vel_pp,accel_p,accel_pp,pos_p,pos_pp
  real(8) :: ljpe
  real(8) :: ke,tot_energy
  real(8) :: pressure_temp
  real(8) :: pressure
  real(8) :: vol
  real(8) :: temp_m

  vel_p    = 0.d0
  vel_pp   = 0.d0
  accel_p  = 0.d0
  accel_pp = 0.d0
  pos_p    = 0.d0
  pos_pp   = 0.d0
  pressure_temp = 0.d0
  ke       = 0.d0
  ljpe     = 0.d0

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
      
      !Calculate Force
      ljpe = 0.d0
      call force_calc(pos_prtl,numprtl,boundry,accel_prtl,ljpe,pressure_temp)   !Calculate acceleration  for t+1
      
      !Update Velocity
      vel_prtl = vel_prtl + 0.5*accel_prtl*dt                !Calculate vel for t+1
      
      
      !Energy 
      ke = .5d0*sum(vel_prtl**2)
      tot_energy = ljpe + ke
      write(51,*),ii,tot_energy,ke,ljpe
      
      !Temperature Writen to file
      temp_m = ke/(1.5d0*dble(numprtl-1))
      write(53,*),ii,temp_m
      
      !Pressure
      vol = boundry**(3.d0)
      pressure = (numprtl*temp_m - (1.d0/(3.d0*temp_m))*pressure_temp)/vol
      write(52,*),ii,pressure

      !Write position to file
      call xyz_file_printout(pos_prtl,numprtl,ii)

      !Scale Velocity
      if (mod(ii,20) == 0 .and. (ii .lt. 200))  then 
          call scale_vel(vel_prtl,numprtl,ke,temp)
          print*,'calling Scale Velocity'
      end if


  end do 
end subroutine!}}}

subroutine scale_vel(vel_prtl,numprtl,ke,temp)!{{{
!http://www.pages.drexel.edu/~cfa22/msim/node33.html
!http://dasher.wustl.edu/bio5476/lectures/lecture-06.pdf  !Input Variable
 
  !Input Variables
  real(8),dimension(numprtl,3) :: vel_prtl
  real(8) :: temp
  real(8) :: ke 
  integer :: numprtl
  
  !internal Variables
  real(8) :: T  !scaled Velocity scalefactor
  real(8) :: alpha
  
  T = ke/(1.5d0*dble(numprtl-1))
  alpha = sqrt(T/temp)
  vel_prtl = vel_prtl/alpha

end subroutine!}}}


