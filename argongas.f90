program ArgonGas

  implicit none
 
! Iniciar variables
! Iniciar matrices
! -posicion
! -momento
! -fuerza

integer :: N = 32, L = 2
real(8) :: T = 1.0,  

REAL, DIMENSION(3,N) :: Posiciones, Momento, Fuerza

! Condiciones iniciales
! -posiciones
! -momento
!  -- rejection

! Iteraciones
! -Calculo de Fuerza
!  -- F = 24*E*[2/r^(13) - 1/r^(7)]
! -Cambiar posiciones
!  -- x = x_i + v*delta_t + 0.5*(F/m)*(delta_t)^2 
!  -- B.C.'s
! -Cambiar momento
!  -- v = F/m*delta_t + v_i

! -Visualization


end program


