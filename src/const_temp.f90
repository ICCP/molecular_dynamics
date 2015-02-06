module const_temp

  private

  public constant_temperature

contains

  subroutine constant_temperature(particles, temp_target, velocity)
    integer, intent(in) :: particles
    real(8), intent(out) :: temp_target
    real(8), intent(inout) :: velocity(3,particles)
    real(8) :: temp_final
    integer :: i

temp_final = 0.0

    do i = 0, particles
       temp_final = temp_final + dot_product(velocity(:,i),velocity(:,i))
    end do
    temp_final = temp_final/(3*(particles-1))
    velocity(:,:) = velocity(:,:)*sqrt(temp_target/temp_final)

  end subroutine constant_temperature

end module const_temp
