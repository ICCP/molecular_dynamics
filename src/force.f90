module force
  
  private
  
  public calculate_force
  public initialize_force

contains

  subroutine initialize_force(forces, N)

    integer, intent(in) :: N
    real(8), intent(inout) :: forces(3,N)

    do i = 1, N
       do j = 1, 3 
          forces(j,i)=0
       end do
    end do
    
  end subroutine initialize_force

  subroutine calculate_force(forces, N, positions, boxes, ener_pot)

    integer, intent(in) :: N, boxes
    real(8), intent(out) :: forces(3,N), ener_pot
    real(8), intent(in) :: positions(3,N)

    real(8) :: distance(3), max, F, rsq
    integer ::i, j

    max = 2.**(2./3)*boxes
    ener_pot =  0

    forces = 0d0

    !write(*,*) 2*max
    !read(*,*)
    do i=1, N - 1
       do j=i+1, N
          distance(:) = positions(:,j)-positions(:,i)       
          !calculates the diference in position
          distance(:) = distance(:) - Nint(distance(:)/(max))*max
          !if this diferences divided by the length of the box (distance > max)
          !then it subtracts L, if is smaller Nint is 0 and doesn't do anything
          rsq = dot_product(distance, distance)
          if(rsq<(3.2**2)) then
             !if (rsq .lt. 1) write(*,*) "====", i, j, rsq
             F =- 24*(2/(rsq**7) - 1/(rsq**4))
             ener_pot = ener_pot + 4*(1/(rsq**6) - 1/(rsq**3))
             !calculates the force
             forces(:,i)=forces(:,i)+distance(:)*F
             forces(:,j)=forces(:,j)-distance(:)*F
          end if
       end do
    end do

  end subroutine calculate_force

end module force
