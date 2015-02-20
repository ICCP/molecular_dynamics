module pressure

  use global

  implicit none

  private
  
  public mean_and_std_dev
  public calc_print_mean_press

contains
  
  subroutine mean_and_std_dev(A, time_sections, mean, std)
    
    real(8), intent(in) :: A(number_timesteps)
    integer, intent(in) :: time_sections
    real(8), intent(out) :: mean
    real(8), intent(out) :: std
    
    real(8) :: mean_array(int((number_timesteps-time_cut)/time_sections))
    real(8) :: mean_prov
    integer :: i,j

    mean_array = 0._8
    mean_prov = 0._8
    i = 1
    mean = 0._8
    std = 0._8
 
    do j = time_cut,  number_timesteps
       mean_prov = mean_prov + A(j)
 !      mean_array(i) = mean_array(i) + A(j)
       if (modulo(j-time_cut+1, time_sections) == 0) then
          mean_prov = mean_prov / time_sections
  !        mean_array(i) = mean_array(i) / time_sections
          mean_array(i) = mean_prov
          mean_prov = 0._8
          i = i+1
       end if
    end do

    mean = sum(mean_array)/(i-1)
    std = dot_product(mean_array, mean_array)/(i-1)
    std = std-(mean)**2
    std = sqrt(std/(i-1))

  end subroutine mean_and_std_dev

  subroutine calc_print_mean_press(mean_pressure, error_pressure)

    real(8), intent(inout) :: mean_pressure, error_pressure
    
    mean_pressure = mean_pressure/(3*num_particles*temp_target)
    mean_pressure = 1 + mean_pressure
    mean_pressure = mean_pressure + 16*PI*num_particles/(3*length**3*temp_target)*(2/3*3.2**(-9)-3.2**(-3))
    error_pressure = error_pressure/(3*num_particles*1.095)
    
    print*, "Average preassure:", mean_pressure
    print*, "Standar error of mean:", error_pressure

  end subroutine calc_print_mean_press

end module pressure
