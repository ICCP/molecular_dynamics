module pressure

  use global

  implicit none

  private
  
  public mean_and_std_dev

contains
  
  subroutine mean_and_std_dev(A, time_sections, mean, std)
    
    real(8), intent(in) :: A(number_timesteps)
    integer, intent(in) :: time_sections
    real(8), intent(out) :: mean
    real(8), intent(out) :: std
    
    real(8) :: mean_array(int((number_timesteps-700)/time_sections))
    integer :: i,j
    
    mean_array = 0._8
    i=1
    mean = 0.0
    std = 0.0
    
    do j = 700,  number_timesteps
       mean_array(i) = mean_array(i) + A(j)
       if (modulo(j-700+1, time_sections) == 0) then
          mean_array(i) = mean_array(i) / time_sections
          i = i+1
       end if
    end do
    
    mean = sum(mean_array)/(i-1)
    std = dot_product(mean_array, mean_array)/(i-1)
    std = std-(mean)**2
    std = sqrt(std/(i-1))
    
  end subroutine mean_and_std_dev

end module pressure
