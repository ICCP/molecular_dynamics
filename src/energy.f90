module energy

  use global

  implicit none

  private

  public write_energy_pressure

contains

  subroutine write_energy_pressure()
    
    integer :: i

    do i = 1, number_timesteps

       write(11,*) i, kin_energy(i)
       write(12,*) i, pot_energy(i)
       write(13,*) i, kin_energy(i)+pot_energy(i)
       write(16,*) i, pressures(i)

    end do

  end subroutine write_energy_pressure

end module energy
