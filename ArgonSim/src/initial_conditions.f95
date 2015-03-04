module initial_conditions

  implicit none
  !private init_random_seed
  public max_boltz, fcc_lattice

contains
  
  subroutine fcc_lattice(sigma, r, cube_side)
    
    real(16), intent(in) :: sigma
    integer(8), intent(in) :: cube_side
    real(16), intent(inout), dimension(:,:) :: r
    ! Define the normal vectors
    integer(8), parameter :: E1(3) = (/ 1, 0, 0 /)
    integer(8), parameter :: E2(3) = (/ 0, 1, 0 /)
    integer(8), parameter :: E3(3) = (/ 0, 0, 1 /)
    real(16), parameter :: CUBE_CORNER(3) = (/ 0, 0, 0 /)
    
    real(16) :: corner_pos(3)
    integer :: i, j, k, n_aux, n
    real(16) :: a

    ! The cube side (a) is given by 2^(2/3) * sigma
    ! in order for the particles on the faces be 
    ! 2^(1/6) * sigma (minimum of the potential)
    ! from the corner particles.
    n = size(r, 2)
    a = 2d0**(2d0/3) * sigma
    print *, "Cube side: ", cube_side
    n_aux = 1
    ! For each cube side we iterate over the corners
    do i = 1, cube_side
      do j = 1, cube_side
        do k = 1, cube_side
          corner_pos = CUBE_CORNER &
           + a * ((i-1)*E1 + (j-1)*E2 + (k-1)*E3)
          r(:, n_aux) = corner_pos
          r(:, n_aux + 1) = corner_pos &
           + 0.5d0 * a * (E1 + E2)
          r(:, n_aux + 2) = corner_pos &
           + 0.5d0 * a * (E1 + E3)
          r(:, n_aux + 3) = corner_pos &
           + 0.5d0 * a * (E2 + E3)
          n_aux = n_aux + 4
        end do
      end do
    end do
    print *, "N_Lattice: ", size(r, 1)
  end subroutine fcc_lattice

  subroutine max_boltz(m, temp, v)
  ! We use Gaussian random to generate the particle momenta 
  ! in each direction for N-1 particles. The last particle 
  ! momenta is such that the sum is zero

    real(16), intent(in) :: temp, m
    real(16), intent(inout), dimension(:,:) :: v
    integer :: i, j, n
    real(16) :: p_sum, p_rnd, var_p

    var_p = temp/m
    N = size(v, 2)

    do i = 1, 3
      p_sum = 0
      do j = 1, N
        call box_muller(p_rnd)
        v(i, j) = var_p*p_rnd
      end do
      !v(i, :) = v(i, :) - sum(v(i, :))/N
      print *, "Mean: ", i, sum(v(i, :))/N
     end do
 	  
  end subroutine max_boltz

  subroutine box_muller(gauss_number)
	! Implementation of polar Box-Muller algorithm
	! to generate Gaussian Random numbers

    real(16) :: x(2)
    real(16), parameter :: TWOPI = 8*atan(1d0)
    real(16), intent(out) :: gauss_number

    call init_random_seed
    call random_number(x)
    !print *, "Random number: ", x
    gauss_number = sqrt(-2*log(x(1))) * cos(TWOPI * x(2))

  end subroutine box_muller

	! Function to initialize the seed of the random
	! number generator

  subroutine init_random_seed()

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

end module
