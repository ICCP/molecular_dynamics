module pressure

  use global

  implicit none

  private
  
  public mean_and_error
  public calc_print_mean_press

contains
  
  subroutine mean_and_error(A, time_sections, mean, error)
    
    !$$ In this subroutine we calculate the mean and 
    !$$ the variance of a correlated array by diving it into boxes
    real(8), intent(in) :: A(number_timesteps)
    integer, intent(in) :: time_sections
    real(8), intent(out) :: mean
    real(8), intent(out) :: error
    
    real(8) :: mean_array(int((number_timesteps+1-time_cut)/time_sections))
    integer :: i,j, final_steps

    mean_array = 0._8
    i = 1
    final_steps = modulo(number_timesteps-time_cut+1, time_sections)

    !$$ Here we go from after the thermostat to the end
    do j = time_cut,  number_timesteps-final_steps
       !$$ We sum all the values in a box together
       mean_array(i) = mean_array(i) + A(j)
       !$$ If we are at the end of a box then:
       if (modulo(j-time_cut+1, time_sections) == 0) then
          !$$ We calculate the mean of that box and move 
          !$$ to the next one
          mean_array(i) = mean_array(i) / time_sections
          i = i+1
       end if
    end do

    !$$ At the end the total mean of the boxes and the standard
    !$$ error of the mean of this boxes are given back
    mean = sum(mean_array)/(i-1)
    error = sum((mean_array-mean)**2)/(i-1)**2
    error = sqrt(error)
    
  end subroutine mean_and_error

  subroutine calc_print_mean_press()

    !$$ In this subroutine we calculate the pressure and 
    !$$ the error of the mean pressure
    real(8) :: mean_pressure, error_pressure
    integer :: time_sec = 300
    
    call mean_and_error(pressures, time_sec, mean_pressure, error_pressure)

    !$$ We calculate the pressure following the virial equation
    mean_pressure = mean_pressure/(3*num_particles*temp_target)
    mean_pressure = 1 + mean_pressure
    mean_pressure = mean_pressure + 16*PI*num_particles/(3*length**3*temp_target)*(2/3*r_cutoff**(-9)-r_cutoff**(-3))

    error_pressure = error_pressure/(3*num_particles*temp_target)
    
    print*, mean_pressure, error_pressure, r_cutoff, "Pressure"

  end subroutine calc_print_mean_press

end module pressure
