module pair_correl

  use global

  implicit none

  private

  public print_pair_corre
  public calc_pair_corre

contains

  subroutine calc_pair_corre(rsq)

    real(8), intent(in) :: rsq

    integer :: corre_dist

    if (flag == 0) then
       corre_dist = nint(sqrt(rsq)*100)
       if(corre_dist <= nint(length/2*100) ) then
          pair_corr(corre_dist) = pair_corr(corre_dist)+1/(4*PI*(0.01**3)*(corre_dist**2))
          ! 0.1 is the delta r
       end if
    end if

  end subroutine calc_pair_corre

  subroutine print_pair_corre()

    integer :: j
 
    pair_corr = pair_corr*2.0*(length**3)
    pair_corr = pair_corr/(num_particles*(num_particles-1))
    pair_corr = pair_corr/(number_timesteps-time_cut)
    
    do j = 1, nint(length*0.5*100)
       write(15,*) j/100.0, pair_corr(j)
    end do
    
  end subroutine print_pair_corre

end module pair_correl
