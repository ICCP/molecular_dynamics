module pressure

  use global

  implicit none

  private
  
  public mean_and_var
  public calc_print_mean_press

contains
  
  subroutine mean_and_var(A, time_sections, mean, variance)
    
    real(8), intent(in) :: A(number_timesteps)
    integer, intent(in) :: time_sections
    real(8), intent(out) :: mean
    real(8), intent(out) :: variance
    
    real(8) :: mean_array(int((number_timesteps+1-time_cut)/time_sections))
    integer :: i,j, final_steps

    mean_array = 0._8
    i = 1
    mean = 0._8
    variance = 0._8
    final_steps = modulo(number_timesteps-time_cut+1, time_sections)

    do j = time_cut,  number_timesteps-final_steps
       mean_array(i) = mean_array(i) + A(j)
       if (modulo(j-time_cut+1, time_sections) == 0) then
          mean_array(i) = mean_array(i) / time_sections
          i = i+1
       end if
    end do

    mean = sum(mean_array)/(i-1)
    variance = dot_product(mean_array, mean_array)/(i-1)
    variance = variance-(mean)**2
    !print*, "siguiente"
    
  end subroutine mean_and_var

  subroutine calc_print_mean_press()

    real(8) :: mean_pressure, error_pressure
    integer :: time_sec = 300
    
    call mean_and_var(pressures, time_sec, mean_pressure, error_pressure)

    mean_pressure = mean_pressure/(3*num_particles*temp_target)
    mean_pressure = 1 + mean_pressure
    mean_pressure = mean_pressure + 16*PI*num_particles/(3*length**3*temp_target)*(2/3*3.2**(-9)-3.2**(-3))
    error_pressure = sqrt(error_pressure*time_sec/(number_timesteps-time_cut))
    error_pressure = error_pressure/(3*num_particles*1.095)
    
    print*, "Average preassure:", mean_pressure
    print*, "Standar error of mean:", error_pressure

  end subroutine calc_print_mean_press

end module pressure
