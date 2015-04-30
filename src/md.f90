program molecular_dynamics
  use iso_fortran_env

  implicit none

  integer :: N_atoms, N_cells, it, tsteps, itrack, ip, jp
  real(8) :: lattice_size, cell_dx, pi, Temp, dt, time
  real(8) :: Energy, dr(3), rdist, force(3), VLJ
  real(8), allocatable :: r(:,:), vel(:,:), acc(:,:)        ! current time step
  real(8), allocatable :: rm(:,:), velm(:,:), accm(:,:)     ! prev time step
  real(8), allocatable :: rm2(:,:), velm2(:,:), accm2(:,:)  ! prev-prev time step

  write(output_unit,*) '====== ICCP Project 2: Molecular dynamics ======'
  pi = 4.d0 * datan(1.d0)

  call init_random_seed
  call read_inputs
  call init_lattice
  call init_velocity

  time = 0.d0

  ! Integrate over time using Velocity Verlet
  do it = 1,tsteps
    time = dble(it) * dt
    rm2 = rm
    velm2 = velm
    accm2 = accm
    rm = r
    velm = vel
    accm = acc

    ! Step [1]: partial update of velocity
    vel = velm + 0.5d0 * dt * accm
    
    ! Step [2]: update position
    r = rm + dt * velm
    r = r + 0.5d0*(dt*dt) * accm

    ! apply periodic BC
    r = mod(r,lattice_size)

    ! Step [3]: update acceleration
    acc = 0.d0
    VLJ = 0.d0 ! Lennard-Jones potential
    do ip = 1,N_atoms-1
      do jp = ip+1,N_atoms
        dr = r(:,ip) - r(:,jp)
        dr = dr - dble(nint(dr/lattice_size))*lattice_size
        rdist = sum(dr * dr)
        write(77,*) rdist
        VLJ = VLJ + 4.d0 * (rdist**(-6.d0) - rdist**(-3.d0))
        force = (48.d0 * rdist**(-7.d0) - 24.d0 * rdist**(-4.d0)) * dr
        acc(:,ip) = acc(:,ip) + force
        acc(:,jp) = acc(:,jp) - force
      enddo
    enddo

    ! Step [4]: finish velocity update
    vel = vel + 0.5d0 * dt * acc

    Energy = VLJ + 0.5d0 * sum(sum(vel * vel,2),1)
    print *,it,'Energy',Energy

    ! update acceleration/force, and update velocity
    if (mod(it,20) == 0) then
      call scale_velocity
    endif

    write(71,*) it, time, r(:,itrack)
    write(72,*) it, time, vel(:,itrack)
    write(73,*) it, time, acc(:,itrack)
    write(74,*) it, time, Energy
  enddo
  write(76,*) Temp,Energy

write(output_unit,*) 'Done'
call exit(0)

contains
!---------------------------------------------------------------------
! Rescale velocity
!---------------------------------------------------------------------
subroutine scale_velocity
  implicit none

  real(8) :: Erg, fac, Tvel

  Erg = 0.5d0 * sum(sum(vel*vel,2),1)
  Tvel = Erg / (1.5d0 * dble(N_atoms-1))

  ! Write temperature to fort.75
  write(75,*) it,time,Tvel
  
  ! Stop rescaling after 500 time steps
  if (it > 500) return

  fac = sqrt(Tvel / Temp)
  vel = vel / fac

  return
end subroutine scale_velocity
!---------------------------------------------------------------------
! Initialize velocity, acceleration arrays
!---------------------------------------------------------------------
subroutine init_velocity
  implicit none

  integer :: i, j
  real(8) :: tmp, tmp2, cmv(3)

  allocate(vel(3,N_atoms))
  allocate(velm(3,N_atoms))
  allocate(velm2(3,N_atoms))
  allocate(acc(3,N_atoms))
  allocate(accm(3,N_atoms))
  allocate(accm2(3,N_atoms))
  vel = 0.d0
  velm = 0.d0
  velm2 = 0.d0
  acc = 0.d0
  accm = 0.d0
  accm2 = 0.d0

  do i = 1,N_atoms
    do j = 1,3
      call random_number(tmp)
      tmp2 = dexp(-tmp*tmp/2.d0) * dsqrt(Temp/2.d0*pi)
      vel(j,i) = tmp2
    enddo
  enddo

  call mass_vel_correction(vel,cmv)
  call scale_velocity

  return
end subroutine init_velocity
!---------------------------------------------------------------------
! Reset the center of mass velocity to zero
!---------------------------------------------------------------------
subroutine mass_vel_correction(v,cmv)
  implicit none

  real(8) :: v(3,N_atoms), cmv(3)
  integer :: i

  cmv = 0.d0

  ! Compute the CMV by averaging
  do i = 1,N_atoms
    cmv = cmv + v(:,i) / dble(N_atoms)
  enddo

  ! Subtract the CMV from the velocities
  do i = 1,N_atoms
    v(:,i) = v(:,i) - cmv
  enddo

  return
end subroutine mass_vel_correction
!---------------------------------------------------------------------
! Initialize the fcc lattice
! Credit to http://www.pa.msu.edu/~duxbury/courses/phy480/fcc.f90
!---------------------------------------------------------------------
subroutine init_lattice
  implicit none

  integer :: i, j, k, l, cnt
  real(8) :: r0(3,4)

  allocate(r(3,N_atoms))
  allocate(rm(3,N_atoms))
  allocate(rm2(3,N_atoms))
  r(:,:) = 0.d0
  rm(:,:) = 0.d0
  rm2(:,:) = 0.d0

  r0(:,:) = cell_dx / 2.d0
  r0(:,1) = 0.d0
  r0(3,2) = 0.d0
  r0(2,3) = 0.d0
  r0(1,4) = 0.d0

  cnt = 0
  do i = 1,N_cells
    do j = 1,N_cells
      do k = 1,N_cells
        do l = 1,4
          cnt = cnt + 1
          r(1,cnt) = r0(1,l) + dble(i-1)*cell_dx
          r(2,cnt) = r0(2,l) + dble(j-1)*cell_dx
          r(3,cnt) = r0(3,l) + dble(k-1)*cell_dx
        enddo
      enddo
    enddo
  enddo

  ! write lattice points to file for plotting
  open(unit=70,file='lattice.out',status='unknown')
  do i = 1,N_atoms
    write(70,'(I10,3ES24.15)') i,r(:,i)
  enddo
  close(70)

  return
end subroutine init_lattice
!---------------------------------------------------------------------
subroutine read_inputs
  implicit none

  open(unit=69,file='inputs.inp',action='read')
  read(69,*) N_cells  ! # cells in one direction
  read(69,*) Temp     ! Temperature
  read(69,*) tsteps, dt ! Simulation duration and time step
  read(69,*) itrack   ! track which particle
  close(69)

  N_atoms = 4 * N_cells**3
  lattice_size = 2.d0 ** (2.d0/3.d0)*dble(N_cells)
  cell_dx = lattice_size / dble(N_cells)
  
  write(output_unit,*) 'Number of cells:',N_cells**3
  write(output_unit,*) 'Number of atoms:',N_atoms
  write(output_unit,*) 'Lattice size:',lattice_size
  write(output_unit,*) 'Cell size:',cell_dx
  write(output_unit,*) 'Temperature:',Temp

  return
end subroutine read_inputs
!---------------------------------------------------------------------
! Stream data from /dev/urandom to use as seed fodder for RNG
! Code taken from:
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!---------------------------------------------------------------------
subroutine init_random_seed
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
  form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
    read(un) seed
    close(un)
  else
  ! Fallback to XOR:ing the current time and pid. The PID is
  ! useful in case one launches multiple instances of the same
  ! program in parallel.
    call system_clock(t)
    if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
      + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
      + dt(3) * 24_int64 * 60 * 60 * 1000 &
      + dt(5) * 60 * 60 * 1000 &
      + dt(6) * 60 * 1000 + dt(7) * 1000 &
      + dt(8)
    end if
    pid = getpid()
    t = ieor(t, int(pid, kind(t)))
    do i = 1, n
      seed(i) = lcg(t)
    end do
  end if
  call random_seed(put=seed)

  return
end subroutine init_random_seed
!---------------------------------------------------------------------
! Note: This routine from the same source as init_random_seed
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
!---------------------------------------------------------------------
function lcg(s)
  integer :: lcg
  integer(int64) :: s
  if (s == 0) then
  s = 104729
  else
  s = mod(s, 4294967296_int64)
  end if
  s = mod(s * 279470273_int64, 4294967291_int64)
  lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
!---------------------------------------------------------------------
end program molecular_dynamics
