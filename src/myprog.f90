program ArgonGas

  implicit none

  integer, parameter :: length = 3, particles = 32
  integer, parameter :: temperature = 300
  real(8) :: position(length,32), velocity(length,32)


  call initialize_position(position, length)
  call initialize_velocity(velocity, length, particles)
      
  write(*,"(3E10.2)") position !3 integers of width 3
  print *, "Velocidades"
  write(*,"(3E10.2)") velocity

contains

  subroutine initialize_position(position, L)

    real(8), intent(inout) :: position(3,32)
    integer, intent(in) :: L
    integer :: i, j, k, n=1

    do i = 1, L-1
       do j = 1, L-1
          do k = 1, L-1
             position(1,n) = i            !First particle in origin
             position(2,n) = j
             position(3,n) = k
             n = n+1
             position(1,n) = i+0.5        !Second particle in face k=1
             position(2,n) = j+0.5
             position(3,n) = k
             n = n+1
             position(1,n) = i             !Third in face i=1
             position(2,n) = j+0.5
             position(3,n) = k+0.5
             n = n+1
             position(1,n) = i+0.5        !Fourth in face j=1
             position(2,n) = j
             position(3,n) = k+0.5
             n = n+1
          end do
       end do
    end do
   
  end subroutine initialize_position
   
  subroutine initialize_velocity(velocity, L, N)

    real(8), intent(inout) :: velocity(3,32)
    integer, intent(in) :: L, N

    integer :: i,j

    do i = 1, L
       do j = 1, N
          velocity(i,j) = gaussRandom()
       end do
    end do

  end subroutine initialize_velocity

  real(8) function gaussRandom() result(rand_vel)


  end function gaussRandom
 

end program
