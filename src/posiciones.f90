module posiciones
  implicit none
  
  private

  public change_posiciones

contains
  
  subroutine change_posiciones(forces, positions, velocity, N, boxes, timestep)

    real(8), intent(in) :: forces(3,32)
    real(8), intent(inout) :: positions(3,32)
    real(8), intent(in) :: velocity(3,32)
    integer, intent(in) :: N, boxes, timestep

    integer :: i, j

    do i = 1, N
       do j = 1, 3
       positions(j,i) = positions(j,i) + velocity(j,i)*timestep + (0.5)*forces(j,i)*(timestep)**2
       end do
    end do 

  end subroutine

end module
