! PHY 480 Project 2
Program lattice
  Implicit None
  real(8):: hold_temp,temp_out,delta_t
  integer max_t, time
  real(8):: mass,potential_energy, kinetic_energy
  real(8):: pressure
  real(8):: length, cell_size
  integer ::cell_dim,N,counter
  real(8), allocatable :: pos(:,:)
  real(8), allocatable :: momenta(:,:)
  real(8), allocatable :: force(:,:)
  

  delta_t = .004d0
  max_t = 5000
  cell_size = 2d0**(2d0/3d0)
  hold_temp = 1d0
  mass = 1d0
  cell_dim = 5
  length = cell_size*cell_dim
  N = 4*cell_dim**3
  allocate( pos(3,N) )
  allocate( momenta(3,N) )
  allocate( force(3,N) )
  force = 0d0
  potential_energy = 0d0
  kinetic_energy = 0d0
  print *, length
  call init_random_seed()
  call position_initializer(N,cell_dim,pos)
  call momentum_initializer(hold_temp,N,mass,momenta)
  call calc_energy(pos,momenta,length,N,potential_energy,kinetic_energy)
  open(unit=20,file = "energy.dat")
  open(unit=19,file = "pressure.dat")
  open(unit=18,file = "momentum.dat")
  
  do time = 1, max_t

     call update(force,pos,momenta,length,N,delta_t,potential_energy,kinetic_energy,hold_temp,temp_out,pressure)
     if (time < 1000 .and. modulo(time,40) == 0) then
        momenta = momenta*sqrt(hold_temp/temp_out)
     end if

     write(20,*) (delta_t*time)," ", (potential_energy+ kinetic_energy), temp_out, kinetic_energy
     write(19,*) delta_t*time, " ",pressure
     write(18,*) delta_t*time," ",sum(momenta(1,:)), " ",sum(momenta(2,:))," ",sum(momenta(3,:))
  end do
  call correlation_new(pos,length,N)
  close(20)     
  close(19)
  close(18)
end program

subroutine correlation_new(position,length,N)
  integer, intent(in) :: N
  real(8), dimension(3,N), intent(in)::position
  real(8), intent(in) :: length
  real(8), dimension(nint(length/2*100)) :: corr_histo
  real(8) :: separation(3),total,new_total,prefactor
  real(8) :: pi
  integer:: counter1,counter2,counter3,r
  total =0d0
  new_total =0d0
  pi = 4*atan(1d0)
  prefactor = 2*(length**3)/(N*(N-1)*4*PI*.01**3)
  print *, length
  print *, prefactor
  do counter1 = 1,N-1
     do counter2 = counter1 + 1, N

        separation = position(:,counter1) - position(:,counter2)
        separation = separation - nint(separation/length)*length
        r = nint(sqrt(dot_product(separation,separation))*100)
        !print *, counter1, counter2, r
        if (r<=nint(length*.5*100)) then
           corr_histo(r) = corr_histo(r) + prefactor/(r**2)
        end if
     end do
  end do

  open(unit=1,file="correlation.dat")
  do counter3 = 1,nint(length*.5*100)
     write(1,*) counter3," ",corr_histo(counter3)
  end do
  

  close(1)
  
end subroutine correlation_new


subroutine calc_energy(position,velocity,length,N,potential,kinetic)
  integer, intent(in) :: N
  real(8), intent(in),dimension(3,N) :: position,velocity
  real(8), intent(in):: length
  real(8), intent(out) :: potential, kinetic
  real(8) :: separation(3), sep_squared
  integer counter1, counter2
  real(8) :: vel_sq(N)
  
  separation = 0d0
  potential = 0d0
  kinetic = 0d0
  do counter1 = 1, N-1
     do counter2 = counter1 + 1, N
        separation = position(:,counter1) - position(:,counter2)
        separation = separation - NINT(separation/length)*length
        sep_squared = dot_product(separation,separation)
        if (sep_squared < (3.2**2)) then
           potential = potential + 4*(1/(sep_squared**6) - 1/(sep_squared**3))
        end if
     end do
  end do
  kinetic = .5*sum(velocity**2)
end subroutine calc_energy

subroutine calc_force(position, force, length, N,temp,pressure)
  real(8), intent(in) ::  length, temp
  integer, intent(in) :: N
  real(8), intent(in),dimension(3,N) :: position
  real(8), intent(out),dimension(3,N):: force
  real(8), intent(out) :: pressure
  real(8) :: separation(3), sep_squared, f(3)
  integer counter1, counter2
  real(8) :: total=0
  separation = 0d0
  force = 0d0
  f = 0d0
  pressure = N*temp/(length**3)
  do counter1 = 1, N-1
     do counter2 = counter1 + 1, N
        separation = position(:,counter1) - position(:,counter2)
        separation = separation - NINT(separation/length)*length 
        sep_squared = DOT_PRODUCT(separation,separation)
        if (sep_squared < (3.2**2)) then
           f = (48/(sep_squared**7)-24/(sep_squared**4))*separation
           force(:,counter1) = force(:,counter1) + f
           force(:,counter2) = force(:,counter2) - f
        end if
        pressure = pressure + sqrt(dot_product(f,f)*sep_squared)/(3d0*temp*length**3)

     end do
  end do  
end subroutine calc_force


subroutine update(force,position,velocity,length,N,delta_t,potential_energy,kinetic_energy,temp_in,temp_out,pressure)
  integer, intent(in) :: N
  real(8), intent(inout),dimension(3,N) :: force, position, velocity
  real(8), intent(in) :: length, delta_t,temp_in
  real(8), intent(out) :: kinetic_energy, temp_out
  real(8), intent(inout) :: potential_energy,pressure
  velocity = velocity + force*delta_t*.5
  position = modulo((position + velocity*delta_t),length)
  call calc_force(position, force, length, N,temp_in,pressure)
  velocity = velocity + force*delta_t*.5
  call calc_energy(position,velocity,length,N,potential_energy,kinetic_energy)
  temp_out = sum(velocity**2)/(3*(N-1))
end subroutine update



subroutine get_gaussian_momentum(mass,stand_dev,x)

  real(8), intent(in) :: stand_dev, mass
  real(8), intent(out) :: x


  real(8) :: random_1,random_y, I, x_integ, dx, pi, prefactor, y, test
  integer :: counter, check
  x = 3.0d0  
  check = 0
  pi = 4d0*atan(1.0d0)
  dx = stand_dev/10d0
  x= 0d0
  prefactor = 1d0/(stand_dev)*sqrt(2d0*pi)
  do while (check == 0)
     call random_number(random_1)
     random_1 = random_1
     x_integ = -10d0*stand_dev
     I = 0d0
     do while (x_integ < 10d0*stand_dev)
        I = I + prefactor*exp(-x_integ*x_integ/(2d0*stand_dev*stand_dev))*dx
        if (I > random_1) then
           x = x_integ*mass
           exit
        else
           x_integ = x_integ + dx
        end if
     end do
     call random_number(random_y)
     random_y = 10d0*random_y
     test = prefactor*exp(-x*x/(2d0*stand_dev*stand_dev))
     if (random_y < test) then
        return
        check = 1
     end if
  end do
end subroutine get_gaussian_momentum




subroutine init_random_seed()
  implicit none

  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  call random_seed(size = n)
  allocate(seed(n))
  open(newunit=un, file="/dev/urandom", access="stream",&
       form="unformatted", action="read", status="old", &
       iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if
     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if
  end if
  call random_seed(put=seed)
end subroutine init_random_seed

subroutine position_initializer(N,dim_cells,position) 
  integer,intent(in) :: N,dim_cells
  real(8),intent(out), dimension(3,N) :: position
  integer :: x_count,y_count,z_count,particle_count
  real(8):: cell_width
  cell_width = 2d0**(2d0/3d0)
  !--------------------------------------------------------------
  !initializing positions of particles in room as FCC lattice
  particle_count =1
  do z_count = 0, dim_cells-1
     do y_count =0,dim_cells-1
        do x_count = 0, dim_cells-1
           position(:,particle_count)=(/ (0d0+x_count)*cell_width,(0d0+y_count)*cell_width,(0d0+z_count)*cell_width /)
           position(:,particle_count+1)=(/ (0.5d0+x_count)*cell_width,(0.5d0 +y_count)*cell_width,(0d0+z_count)*cell_width /)
           position(:,particle_count+2)=(/ (0.5d0+x_count)*cell_width,(0d0+y_count)*cell_width,(0.5d0+z_count)*cell_width /)
           if (particle_count+3 > N) then
              print *, "Out of bounds in position array"
              exit
           end if
           position(:,particle_count+3) = (/ (0d0+x_count)*cell_width,(0.5d0+y_count)*cell_width,(0.5d0+z_count)*cell_width /)
           particle_count = particle_count + 4
        end do
     end do
  end do
  !---------------------------------------------------------------
  return
end subroutine position_initializer

subroutine momentum_initializer(temp, N, mass,momentum)
  real(8) :: temp, mass, stand_dev, kb
  integer :: N, counter
  real(8), intent(out), dimension(3,N) :: momentum
  real(8), dimension(3) :: total_momentum
  kb = 1d0
  stand_dev = sqrt(kb*mass*temp)
  !---------------------------------------------------------------
  !initializing velocities of particles according to a Gaussian distribution
  do counter =1,N
     call get_gaussian_momentum(mass,stand_dev,momentum(1,counter))
     call get_gaussian_momentum(mass,stand_dev,momentum(2,counter))
     call get_gaussian_momentum(mass,stand_dev,momentum(3,counter))
  end do
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !finds total momentum in each direction and reduces to get closer to 0
  total_momentum = (/sum(momentum(1,:)), sum(momentum(2,:)), sum(momentum(3,:))/)

  momentum(1,:) = momentum(1,:) - total_momentum(1)/N
  momentum(2,:) = momentum(2,:) - total_momentum(2)/N
  momentum(3,:) = momentum(3,:) - total_momentum(3)/N
  total_momentum = (/sum(momentum(1,:)), sum(momentum(2,:)), sum(momentum(3,:))/)
  print *,total_momentum(1),total_momentum(2),total_momentum(3)
  !---------------------------------------------------------------
  return
end subroutine momentum_initializer
