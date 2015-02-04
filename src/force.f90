module force
  
  private
  
  public calculate_force

contains

  subroutine calculate_force(forces, N, positions, boxes)

    real(8), intent(inout) :: forces(3,32)
    real(8), intent(in) :: positions(3,32)
    integer, intent(in) :: N, boxes

    real(8) :: distance(3), max =2.0**(2.0/3)
    integer ::i, j

    do i = 1, 32
       do j = 1, 3 
          forces(j,i)=0
       end do
    end do

    i=1!, N
       do j=1, N
          distance = sqrt((positions(1,i)-positions(1,j))**2+(positions(2,i)-positions(2,j))**2+(positions(3,i)-positions(3,j))**2)
          !print *, j, distance
          distance = mod(distance, max)
          !print *, distance

       end do
    !end do

  end subroutine calculate_force



end module force
