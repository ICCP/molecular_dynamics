module energy

  use global
  use pressure

  implicit none

  private

  public write_files
  public calc_print_heat_capac

contains

  subroutine write_files()
    
    integer :: i

    do i = time_cut, number_timesteps

       write(11,*) i, vel_sqr(i)/(2*num_particles)
       write(12,*) i, pot_energy(i)/(num_particles)
       write(13,*) i, (vel_sqr(i)/2+pot_energy(i))/num_particles
       write(14,*) i, vel_sqr(i)/(3*(num_particles-1))
       write(16,*) i, pressures(i)

    end do

  end subroutine write_files

  subroutine calc_print_heat_capac()

    real(8) :: mean_kinen, variance_kinen
    integer :: time_sec = 200

    call mean_and_std_dev(pot_energy, time_sec, mean_kinen, variance_kinen)

    variance_kinen = variance_kinen**2*int((number_timesteps-time_cut)/time_sec)

    print*, mean_kinen, variance_kinen
    print*, "Heat Capacitance (cv):", (2/(3*num_particles)-variance_kinen/mean_kinen**2)**(-1)/num_particles

  end subroutine calc_print_heat_capac

end module energy
