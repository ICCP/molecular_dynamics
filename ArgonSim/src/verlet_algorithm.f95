module verlet

implicit none
public verlet_algorithm, total_force, calc_energy

contains
    subroutine verlet_algorithm(EPS, SIGMA, MASS, l, dt, r, v, a, Ep, vsquared, olda, temp_a, pressure, currTemp, kB)
      real(16), intent(in) :: SIGMA, MASS, EPS, dt, kB, l
      real(16), intent(inout), dimension(:,:) :: r, v, a, olda
      real(16), intent(inout) :: vsquared, currTemp, temp_a(3)
      real(16), intent(out) :: Ep, pressure
      real(16) :: temp_r
      integer :: N, i, j
      N = size(r,2)

      call total_force(EPS, SIGMA, MASS, l, r, a, Ep, kB, currTemp, pressure, temp_a)
      do i=1, N
        do j=1, 3
          temp_r = r(j,i)+v(j,i)*dt+(0.5*a(j,i)*MASS)*dt**2d0
          r(j,i) = modulo(temp_r, l)
        end do
      end do

      olda = a
      call total_force(EPS, SIGMA, MASS, l, r, a, Ep, kB, currTemp, pressure, temp_a)

      v = v + 0.5*dt*MASS*(a+olda)
      do i=1,3
        v(i,:) = v(i,:) - sum(v(i,:))/N
      end do
      vsquared = sum(v**2)/N
      !print *, "a",a(:,1),"v",v(:,1),"r",r(:,1)

    end subroutine verlet_algorithm

    subroutine total_force(EPS, SIGMA, MASS, l, r, a, Ep, kB, currTemp, pressure, temp_a) ! n_aux is the particle we're calculating the force on

        real(16), intent(inout), dimension(:,:) :: r, a
        real(16), intent(in) :: SIGMA, MASS, EPS, l
        real(16), intent(inout):: Ep, pressure, currtemp
        real(16) :: pressure_expect, kB
        integer :: N, i, j, k
        real(16) :: diff_aux(3), dist, fij(3), temp_a(3)
        N = size(r,2)
        a(:,:) = 0
        pressure_expect = 0 
        do i = 1, N-1
            do j = (i + 1), N
                diff_aux = r(:, j) - r(:, i)
                diff_aux = diff_aux - nint(diff_aux/l)*l
                dist = dot_product(diff_aux, diff_aux)
                ! No need for force array; just use acceleration / mass (same dimensions, too)
                if (dist < 3.2**2) then
                  temp_a = 4d0*EPS*(12*(SIGMA**12)/(dist**7)-6*(SIGMA**6)/(dist**4))*(diff_aux(:))/MASS
                else
                  temp_a = 0
                end if
                a(:, i) = a(:, i) - temp_a
                a(:, j) = a(:, j) + temp_a
                pressure_expect = pressure_expect + dist*sqrt(dot_product(temp_a, temp_a))
            end do
        end do
        pressure = (N * kB * currTemp + pressure_expect/(3.0*MASS*N)) / l**3.0

    end subroutine total_force

    subroutine calc_energy(v, MASS, Ek, Ep, mom, SIGMA, EPS, r, l)
      real(16), intent(in), dimension(:,:) :: v, r
      real(16), intent(in) :: MASS, l, SIGMA, EPS
      real(16), intent(out) :: Ek, Ep, mom(3)
      real(16) :: diff_aux(3), dist, EnAvg
      integer :: N, i, j

      N = size(v, 2)
      Ep = 0
      Ek = 0
      do i = 1, N-1
        Ek = Ek + 0.5d0*MASS*dot_product(v(:, i), v(:, i))
        do j = (i + 1), N
            diff_aux = r(:, j) - r(:, i)
            diff_aux = diff_aux - nint(diff_aux/l)*l
            dist = dot_product(diff_aux, diff_aux)
            Ep = Ep + 4d0*EPS*((SIGMA**2d0/dist)**(6d0) - (SIGMA**2d0/dist)**(3d0))
            end do
      end do
      Ek = Ek + 0.5d0*MASS*dot_product(v(:, N), v(:, N))
    end subroutine calc_energy

    subroutine pair_correlation(dr, r, r_histo, l)
      real(16), intent(in) :: dr, l
      real(16), intent(in), dimension(:,:) :: r
      real(16), intent(inout), dimension(:) :: r_histo
      real(16) :: diff_aux(3), dist
      integer :: i, N, j, bin, bin_count

      N = size(r,2)
      bin_count = size(r_histo)
      do i = 1, N
        do j = (i + 1), N
          diff_aux = r(:, j) - r(:, i)
          diff_aux = diff_aux - nint(diff_aux/l)*l
          dist = sqrt(dot_product(diff_aux, diff_aux))
          print *, "Dist: ", dist
          print *, "dr: ", dr
          bin = int(dist/dr)
          if (bin > bin_count) then
            print *, "WARNING: Bin > bin_count"
          end if
          print*, "bin: ", bin
          r_histo(bin) = r_histo(bin) + 1
        end do
      end do

      do i = 1, bin_count
        if (r_histo(i) /= 0) then
          r_histo(i) = 3*r_histo(i)/(4*3.1415*((dr)**3)*(((1+i)**3)-i**3))
        end if
      end do
    end subroutine pair_correlation
end module