module const_temp

  implicit none

  private

  public constant_temperature

contains

  subroutine constant_temperature(particles, temp_target, velocity, temp_final)
    integer, intent(in) :: particles
    real(8), intent(in) :: temp_target
    real(8), intent(inout) :: velocity(3,particles)
    real(8), intent(out) :: temp_final
    
    temp_final = sum(velocity**2)
    temp_final = temp_final/(3*(particles-1))
    velocity(:,:) = velocity(:,:)*sqrt(temp_target/temp_final)
  end subroutine constant_temperature
  
end module const_temp
