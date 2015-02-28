module force

  use global
  use pair_correl

  implicit none
  
  private 

  public interac_part

contains

  subroutine interac_part(time)

    !$$ In this subroutine we calculate the interaction between
    !$$ all the particles. We save the forces, the potential
    !$$ energy and the velocity square for further calculations. 
    integer, intent(in) :: time
 
    real(8) :: F, rsq, distance(3)
    integer :: i, j

    forces = 0._8

    do i=1, num_particles - 1
       do j=i+1, num_particles

          !$$ We calulate the diference in position of every particle and
          !$$ apply the boundary conditions to it. 
          distance = position(:,j)-position(:,i)
          distance = distance - Nint(distance/length)*length
          rsq = dot_product(distance, distance)
          !$$ In this moment we calculate the pair correleation function
          call calc_pair_corre(rsq)
          !$$ And if the distance is smaller than the r_cutoff we
          !$$ calculate the force, the potential energy and the pressure
          if(rsq<(r_cutoff**2)) then
             F = 24*(2/(rsq**7) - 1/(rsq**4))
             forces(:,i)=forces(:,i)-distance*F
             forces(:,j)=forces(:,j)+distance*F
             if (time >= time_cut) then          
                pot_energy(time) = pot_energy(time) + 4*(1/(rsq**6) - 1/(rsq**3)) 
                pressures(time) = pressures(time) + rsq*F
             end if
          end if
       end do
    end do

  end subroutine interac_part

end module force
