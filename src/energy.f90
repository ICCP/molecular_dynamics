module energy

  use global
  use pressure

  implicit none

  private calc_cv_coeff

  public calc_energies
  public calc_print_heat_capac

contains

  subroutine calc_energies()
    
    !$$ In this subroutine we calculate the temperature
    !$$ and the potential energy means and errors
    real(8) :: mean_temp, error_temp
    real(8) :: mean_poten, error_poten
    integer, parameter :: time_sec = 300

    call mean_and_error(vel_sqr, time_sec, mean_temp, error_temp)
    call mean_and_error(pot_energy, time_sec, mean_poten, error_poten)

    mean_temp = mean_temp/(3*(num_particles-1))
    error_temp = error_temp/(3*(num_particles-1))
    print*, mean_temp, error_temp, r_cutoff, "Temperature"

    mean_poten = mean_poten/num_particles
    error_poten = error_poten/num_particles
    print*, mean_poten, error_poten, r_cutoff, "Potential Energy"

  end subroutine calc_energies

  subroutine calc_print_heat_capac()

    !$$ To calculate the specific heat we use the Lebowitz formula
    real(8) :: cv, error_cv
    integer, parameter :: time_sec = 300
    real(8) :: cv_coefficient(number_timesteps)

    !$$ In thgis equation we have a coefficient of the
    !$$ variance divided by the mean. In this function we
    !$$ calculate this coefficient every timestep
    call calc_cv_coeff(cv_coefficient, time_sec)
   
    !$$ Then we calculate the mean and the error of them to 
    !$$ get the mean and error
    call mean_and_error(cv_coefficient, time_sec, cv, error_cv)
    cv = (2.0/(3.0)-cv*num_particles)**(-1)
    error_cv = error_cv*num_particles

    print*, cv, error_cv, r_cutoff, "heat capacity"

  end subroutine calc_print_heat_capac

  subroutine calc_cv_coeff(cv_coeff, time_sec)

    real(8), intent(out) :: cv_coeff(number_timesteps)
    integer, intent(in) :: time_sec

    real(8) :: mean_array(int((number_timesteps+1-time_cut)/time_sec))
    integer :: i, j, final_steps

    !$$ The coefficients are first the kinetic energy per particle
    cv_coeff = vel_sqr/(2*num_particles)

    mean_array = 0._8
    i = 1
    final_steps = modulo(number_timesteps-time_cut+1, time_sec)

    !$$ Then we calculate the average of them in boxes
    do j = time_cut,  number_timesteps-final_steps
       mean_array(i) = mean_array(i) + cv_coeff(j)
       if (modulo(j-time_cut+1, time_sec) == 0) then
          mean_array(i) = mean_array(i) / time_sec
          i = i+1
       end if
    end do

    !$$ And then we calculate the coefficient with the mean of its box
    i=1
    do j = time_cut,  number_timesteps-final_steps
       cv_coeff(j) = (cv_coeff(j)-mean_array(i))/mean_array(i)
       cv_coeff(j) = cv_coeff(j)**2
       if (modulo(j-time_cut+1, time_sec) == 0) then
          i = i+1
       end if
    end do


  end subroutine calc_cv_coeff

end module energy
