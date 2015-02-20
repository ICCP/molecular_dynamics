module pair_correl

  use global

  implicit none

  private

  public print_pair_corre
  public calc_pair_corre

contains

  subroutine calc_pair_corre(flag, rsq, length, pair_corre)

    integer, intent(in) :: flag
    real(8), intent(in) :: rsq, length
    real(8), intent(inout) :: pair_corre(nint(length/2*100))

    integer :: corre_dist

    if (flag == 0) then
       corre_dist = nint(sqrt(rsq)*100)
       if(corre_dist <= nint(length/2*100) ) then
          pair_corre(corre_dist) = pair_corre(corre_dist)+1/(4*PI*(0.01**3)*(corre_dist**2))
       end if
    end if

  end subroutine calc_pair_corre

  subroutine print_pair_corre(pair_corre)

    real(8), intent(inout) :: pair_corre(nint(length/2*100))

    integer :: j
    real(8) :: length
    length = nint((num_particles/4)**(1.0/3))*(4.0/density)**(1.0/3)

    pair_corre = pair_corre*2.0*(length**3)
    pair_corre = pair_corre/(num_particles*(num_particles-1))
    pair_corre = pair_corre/(number_timesteps-700)
    
    do j = 1, nint(length*0.5*100)
       write(15,*) j/100.0, pair_corre(j)
    end do
    
  end subroutine print_pair_corre

end module pair_correl
