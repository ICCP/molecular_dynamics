module pair_correl

  use global

  implicit none

  private

  public print_pair_corre
  public calc_pair_corre

contains

  subroutine calc_pair_corre(rsq)

    !$$ In this subroutine we calculate the pair correlation
    !$$ function by counting how many particles are in a
    !$$ certain distance
    real(8), intent(in) :: rsq

    integer :: corre_dist

    !$$ If the thermostat if of then we calculate it
    if (flag == 0) then
       corre_dist = nint(sqrt(rsq)/delta_r)
       !$$ We assure that the correlation distance is smaller than
       !$$ half of the length
       if(corre_dist <= nint(length/(2*delta_r)) ) then
          !$$ we add up all the particles. We add them already normalized.
          pair_corr(corre_dist) = pair_corr(corre_dist)+1/(4*PI*(delta_r**3)*(corre_dist**2))
       end if
    end if

  end subroutine calc_pair_corre

  subroutine print_pair_corre()

    !$$ Here we print the pair correlation function in a file
    integer :: j
 
    pair_corr = pair_corr*2.0*(length**3)
    pair_corr = pair_corr/(num_particles*(num_particles-1))
    pair_corr = pair_corr/(number_timesteps-time_cut)
    
    do j = 1, nint(length/(2*delta_r))
       write(15,*) j*delta_r, pair_corr(j)
    end do
    
  end subroutine print_pair_corre

end module pair_correl
