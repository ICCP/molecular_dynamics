module time_evol

  use global
  use force
  use energy
  use pressure
  use pair_correl

  implicit none

  private Verlert_alg
  private thermostat
!!$  private plot_particles

  public Time_evolution

contains

  subroutine Time_evolution()

    !$$ Here we do the time evolution of the system
    integer :: i

    !$$ At the begging we need to calculate the first forces
    !$$ that the particles feel
    call interac_part(1)

    do i = 1, number_timesteps
       !$$ For all times we calculate the Verlet Algorithm
       call Verlert_alg(i)
       !$$ We store the information of the square of the velocities
       vel_sqr(i) = sum(velocity**2)
       !$$ And call the thermostat. This will only work after
       !$$ passing time_cut steps
       call thermostat(i)
!!$       call plot_particles()
    end do
    
  end subroutine Time_evolution

  subroutine Verlert_alg(i)

    integer, intent(in) :: i

    !$$ The Verlet Algorithm consist of moving the velocities
    !$$ for have the time_step
    velocity = velocity + forces*time_step/2.0
    !$$ Then update the position with full time_steps
    !$$ and boundary conditions
    position = modulo(position+velocity*time_step, length)
    !$$ Later update the forces with this new positions
    call interac_part(i)
    !$$ And finally do the other half time_step for the velocities
    velocity = velocity + forces*time_step/2.0
    
  end subroutine Verlert_alg

  subroutine thermostat(steps)

    !$$ In this subroutine we apply the thermostat to the system
    integer, intent(in) :: steps

    !$$ If we before the tom_cut we enter the thermostat
    if(steps < time_cut) then
       !$$ Every 20 time_steps we apply the thermostat
       if(modulo(steps,20) == 0) then    
          !$$ Applying a thermostat is only renormalizing the velocity
          !$$ in order to have constant total energy and the system
          !$$ relaxes to the users temperature-parameter
          velocity = velocity*sqrt((3*(num_particles-1)*temp_target/vel_sqr(steps)))
       end if
    !$$ If we just finish the thermostat we take our flag down
    !$$ in order to begin the calculation of the physical quantities
    else if(steps == time_cut) then
       flag = 0
    end if
    
  end subroutine thermostat

  subroutine plot_particles()

    integer :: i

    do i= 1, num_particles
       call setpoint(position(1,i), position(2,i))
    end do

  end subroutine plot_particles
  
end module time_evol
