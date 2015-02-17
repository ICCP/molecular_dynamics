module force

  implicit none
  
  private
  
  public calculate_force

contains

  subroutine calculate_force(forces, N, positions, length, ener_pot, pair_corre, flag, pressure)

    integer, intent(in) :: N, flag
    real(8), intent(out) :: forces(3,N), ener_pot
    real(8), intent(in) :: positions(3,N)
    real(8), intent(in) :: length
    real(8), intent(out) :: pair_corre(nint(length/2*100))
    real(8), intent(out) :: pressure

    real(8), parameter :: PI = 4*atan(1.0)
    real(8) :: distance(3), F, rsq
    integer :: corre_dist
    integer ::i, j

    ener_pot = 0.0
    forces = 0d0
    pressure = 0d0

    do i=1, N - 1
       do j=i+1, N
          distance = positions(:,j)-positions(:,i)       
          !calculates the diference in position
          distance = distance - Nint(distance/(length))*length
          !if this diferences divided by the length of the box (distance > max)
          !then it subtracts L, if is smaller Nint is 0 and doesn't do anything
          rsq = dot_product(distance, distance)

          if (flag == 0) then
             corre_dist = nint(sqrt(rsq)*100)
             if(corre_dist <= nint(length/2*100) ) then
                pair_corre(corre_dist) = pair_corre(corre_dist)+1/(4*PI*(0.01**3)*(corre_dist**2))
             end if
          end if

          if(rsq<(3.2**2)) then
             F =- 24*(2/(rsq**7) - 1/(rsq**4))
             ener_pot = ener_pot + 4*(1/(rsq**6) - 1/(rsq**3))
             forces(:,i)=forces(:,i)+distance*F
             forces(:,j)=forces(:,j)-distance*F
             
             pressure = pressure + rsq*F
          end if


          
       end do
    end do

  end subroutine calculate_force

end module force
