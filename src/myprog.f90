program ArgonGas

  implicit none

  integer, parameter :: length = 3, particles = 32
  integer, parameter :: temperature = 300
  real(8) :: position(3,32), velocity(3,32)


  call initialize_position(position, length)
 ! call initialize_velocity(velocity, length, particles)
      
  write(*,"(3E10.2)") position !3 reals of width 10 and 2 decimals
  print *, "Velocidades" 
  write(*,"(3E10.2)") velocity

contains

  subroutine initialize_position(position, L)

    integer, intent(in) :: L
    real(8), intent(inout) :: position(3,32)
    integer :: i, j, k, n=1

    do i = 0, L-2
       do j = 0, L-2
          do k = 0, L-2
             position(1,n) = i*sqrt(2.0)            !First particle in origin
             position(2,n) = j*sqrt(2.0)
             position(3,n) = k*sqrt(2.0)
             n = n+1
             position(1,n) = i*sqrt(2.0)+1/sqrt(2.0)!Second particle in face k=1
             position(2,n) = j*sqrt(2.0)+1/sqrt(2.0)
             position(3,n) = k*sqrt(2.0)
             n = n+1
             position(1,n) = i*sqrt(2.0)            !Third in face i=1
             position(2,n) = j*sqrt(2.0)+1/sqrt(2.0)
             position(3,n) = k*sqrt(2.0)+1/sqrt(2.0)
             n = n+1
             position(1,n) = i*sqrt(2.0)+1/sqrt(2.0)!Fourth in face j=1
             position(2,n) = j*sqrt(2.0)
             position(3,n) = k*sqrt(2.0)+1/sqrt(2.0)
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
          !velocity(i,j) = gaussRandom()
       end do
    end do

  end subroutine initialize_velocity

 ! real(8) function gaussRandom() result(rand_vel)


 ! end function gaussRandom
 

end program
