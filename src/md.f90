program molecular_dynamics
  use iso_fortran_env

  implicit none

  integer :: N_atoms, N_cells
  real(8) :: lattice_size, cell_dx
  real(8), allocatable :: r(:,:), vel(:,:), acc(:,:)

  write(output_unit,*) '====== ICCP Project 2: Molecular dynamics ======'
  call read_inputs
  call init_lattice

  call exit(0)

contains
!---------------------------------------------------------------------
! Initialize velocity
!---------------------------------------------------------------------
subroutine init_velocity
end subroutine init_velocity
!---------------------------------------------------------------------
! Initialize the fcc lattice
! Credit to http://www.pa.msu.edu/~duxbury/courses/phy480/fcc.f90
!---------------------------------------------------------------------
subroutine init_lattice
  implicit none

  integer :: i, j, k, l, cnt
  real(8) :: r0(3,4)

  allocate(r(3,N_atoms))
  r(:,:) = 0.d0

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
  read(69,*) N_cells ! # cells in one direction
  close(69)

  N_atoms = 4 * N_cells**3
  lattice_size = 2.d0 ** (2.d0/3.d0)*dble(N_cells)
  cell_dx = lattice_size / dble(N_cells)
  
  write(output_unit,*) 'Number of cells:',N_cells**3
  write(output_unit,*) 'Number of atoms:',N_atoms
  write(output_unit,*) 'Lattice size:',lattice_size
  write(output_unit,*) 'Cell size:',cell_dx

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
