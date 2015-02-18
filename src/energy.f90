module energy

  implicit none

  private

  public write_energy

contains

  subroutine write_energy(N, velo, ener_pot, i)

    integer, intent(in) :: N, i
    real(8), intent(in) :: velo(3,N)
    real(8), intent(in) :: ener_pot

    real(8) :: ener_kin

    ener_kin = sum(velo**2*0.5)
    write(11,*) i, ener_kin
    write(12,*) i, ener_pot
    write(13,*) i, ener_kin+ener_pot

  end subroutine write_energy

end module energy
