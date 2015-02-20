module force

  use global
  use pair_correl

  implicit none
  
  private 

  public relation_part

contains

  subroutine relation_part(time)

    integer, intent(in) :: time
 
    real(8) :: F, rsq, distance(3)
    integer :: i, j

    forces = 0._8

    do i=1, num_particles - 1
       do j=i+1, num_particles

          distance = position(:,j)-position(:,i)       
          !calculates the diference in position
          distance = distance - Nint(distance/length)*length
          !if this diferences divided by the length of the box (distance > max)
          !then it subtracts L, if is smaller Nint is 0 and doesn't do anything
          rsq = dot_product(distance, distance)
          call calc_pair_corre(rsq)          
          if(rsq<(3.2**2)) then
             F =- 24*(2/(rsq**7) - 1/(rsq**4))
             forces(:,i)=forces(:,i)+distance*F
             forces(:,j)=forces(:,j)-distance*F
             
             pot_energy(time) = pot_energy(time) + 4*(1/(rsq**6) - 1/(rsq**3))            
             pressures(time) = pressures(time) - rsq*F
          end if

       end do
    end do

  end subroutine relation_part

end module force
