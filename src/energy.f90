module energy

  private

  public calculate_kin_energy

contains

  subroutine calculate_kin_energy(velocity, particles, ener_kin)

    integer, intent(in) :: particles
    real(8), intent(inout) :: velocity(3,particles), ener_kin

    integer :: i

    ener_kin = 0

    do i = 1, particles
       ener_kin = ener_kin +(0.5)*norm2(velocity(:,i))**2
    end do

  end subroutine calculate_kin_energy

end module energy
