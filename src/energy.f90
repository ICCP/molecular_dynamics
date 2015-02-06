module energy

  implicit none

  private

  public calculate_kin_energy

contains

  subroutine calculate_kin_energy(velocity, particles, ener_kin)

    integer, intent(in) :: particles
    real(8), intent(in) :: velocity(3,particles)
    real(8), intent(out) :: ener_kin

    ener_kin = sum(velocity**2*0.5)

!!$    ener_kin = 0
!!$
!!$    do i = 1, particles
!!$       ener_kin = ener_kin +(0.5)*norm2(velocity(:,i))**2
!!$    end do

  end subroutine calculate_kin_energy

end module energy
