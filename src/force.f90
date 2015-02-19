module force

  use global
  use pair_correl

  implicit none
  
  private

  public relation_part

contains

  subroutine relation_part(ener_pot, pair_corre, flag, pressure)

    integer, intent(in) :: flag
    real(8), intent(out) :: ener_pot
    real(8), intent(out) :: pair_corre(nint(length/2*100))
    real(8), intent(out) :: pressure

    real(8) :: distance(3), F, rsq
    integer ::i, j

    ener_pot = 0._8
    forces = 0._8
    pressure = 0._8

    do i=1, num_particles - 1
       do j=i+1, num_particles

          distance = position(:,j)-position(:,i)       
          !calculates the diference in position
          distance = distance - Nint(distance/(length))*length
          !if this diferences divided by the length of the box (distance > max)
          !then it subtracts L, if is smaller Nint is 0 and doesn't do anything
          rsq = dot_product(distance, distance)

          call calc_pair_corre(flag, rsq, length, pair_corre) 

          if(rsq<(3.2**2)) then
             F =- 24*(2/(rsq**7) - 1/(rsq**4))
             forces(:,i)=forces(:,i)+distance*F
             forces(:,j)=forces(:,j)-distance*F
             
             ener_pot = ener_pot + 4*(1/(rsq**6) - 1/(rsq**3))            
             pressure = pressure - rsq*F
          end if
          
       end do
    end do

  end subroutine relation_part

end module force
