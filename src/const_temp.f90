module const_temp

  private

  public constant_temperature

contains

  subroutine constant_temperature(particles, temp_target, velocity, temp_final)
    integer, intent(in) :: particles
    real(8), intent(in) :: temp_target
    real(8), intent(inout) :: velocity(3,particles)
    real(8), intent(inout) :: temp_final
    integer :: i

    do i = 0, particles
       temp_final = temp_final + dot_product(velocity(:,i),velocity(:,i))
    end do
    temp_final = temp_final/(3*particles)
    velocity(:,:) = velocity(:,:)*sqrt(temp_target/temp_final)

  end subroutine constant_temperature

end module const_temp
