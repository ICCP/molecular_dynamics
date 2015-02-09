module force

  implicit none
  
  private
  
  public calculate_force

contains

  subroutine calculate_force(forces, N, positions, length, ener_pot, pair_corre, flag, virial_function)

    integer, intent(in) :: N, flag
    real(8), intent(out) :: forces(3,N), ener_pot
    real(8), intent(in) :: positions(3,N)
    real(8), intent(in) :: length
    real(8), intent(out) :: pair_corre(nint(length*sqrt(3.0)*0.5/0.01))
    real(8), intent(out) :: virial_function

    real(8), parameter :: PI = 4*atan(1.0)
    real(8) :: distance(3), F, rsq
    integer :: corre_dist
    integer ::i, j

    ener_pot = 0.0
    forces = 0d0
    virial_function = 0d0
    do i=1, N - 1
       do j=i+1, N
          distance(:) = positions(:,j)-positions(:,i)       
          !calculates the diference in position
          distance(:) = distance(:) - Nint(distance(:)/(length))*length
          !if this diferences divided by the length of the box (distance > max)
          !then it subtracts L, if is smaller Nint is 0 and doesn't do anything
          rsq = dot_product(distance, distance)

          if (flag == 0) then
             corre_dist = nint(sqrt(rsq)*100)
             pair_corre(corre_dist) = pair_corre(corre_dist)+1/(4*PI*(0.01**3)*(corre_dist**2))
          end if

!!$          if(rsq<(3.2**2)) then
             !if (rsq .lt. 1) write(*,*) "====", i, j, rsq
             F =- 24*(2/(rsq**7) - 1/(rsq**4))
             ener_pot = ener_pot + 4*(1/(rsq**6) - 1/(rsq**3))
             !calculates the force
             forces(:,i)=forces(:,i)+distance*F
             forces(:,j)=forces(:,j)-distance*F
!!$          end if
             virial_function = virial_function + dot_product(distance(:),(distance(:)*F))
       end do
    end do

  end subroutine calculate_force

end module force
