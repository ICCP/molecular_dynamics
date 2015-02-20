module energy

  use global

  implicit none

  private

  public write_files

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

end module energy
