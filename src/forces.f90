    subroutine update_forces(positions, forces, correlation, potential, N, rc, L, correlation_steps, pressure, kB, T, V)
    implicit none

    integer, intent(in) :: N, correlation_steps	! Respectfully the number of particles in the box and number of individual shells in the space 						  	  discretization of the correlation function calculation
    real(8), intent(in) :: positions(N,3), rc, kB, T, V	! The vector containing the position vector of every particle, the critical distance, the temperature and the volume of the box
    real(8), intent(inout) :: forces(N,3), potential, correlation(correlation_steps), pressure
						! The vector containing the force vectors for every particle, the total potential, he 							correlation function (before renormalization, this is done in the plot.correlation_function) and the pressure

!f2py intent(in,out) :: forces, potential, correlation, pressure
    

    real(8) :: dr(3), L, force, a, d, term1, term2
    integer :: i,j,k
    forces = 0._8
    potential = 0._8
    pressure = 0._8
    term1 = 0._8
    a = correlation_steps * 2.0 / L  	! Constant used to increment the correct component in correlation vector

    do i= 1, N 
      do j = i+1, N			
        dr = (positions(i,:)-positions(j,:))	! Distance vector between the particles
        dr = dr - nint(dr/L)*L			! Correction due to periodic boundary conditions
        d = sqrt(sum(dr*dr))			! Calculation of the distance 
        if (d<L/2) then   			! We only go further if the distance is useful for the correlation function calculation (L/2)
		k = nint(a * d)						! Index of the component corresponding to this distance
		correlation(k) = correlation(k) + 1			! We increment the component corresponding to this distance
		if (d<rc) then			! We only go further if the distance is smaller to the critical distance rc 
			force = 24.0*(2.0 / d**14 - 1 / d**8)			! The norm of the force between particle i and j
                	forces(i,:) = forces(i,:) + force*dr 			! The force of particle j acting on i is added to the force vector
			forces(j,:) = forces(j,:) - force*dr  			! The force of particle i acting on j is added to the force vector
			potential = potential + 4 *  (1 / d ** 12 - 1/ d ** 6) 	! The interaction potential is added to the total potential
			term1 = term1 - force*d**2	! Contribution to the virial term for pressure calculation
		end if	
        end if
      end do
    end do
    term2 = (48/9)/rc**9 - (24/3)/rc**3    	! Pair-correlation function contribution to the pressure
    pressure = kb * T * N / V - term1/(3*V*N) - term2*2*3.14*N**2/(3*V**2)
    end subroutine

!To compile and link to Python code: 
!gfortran -c forces.f90
!f2py -c --f90flags='-march=native -ffast-math' forces.f90 -m f90force
