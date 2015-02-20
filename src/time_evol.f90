module time_evol

  use global
  use force
  use energy
  use pressure
  use pair_correl

  implicit none

  private Verlert_alg
  private thermostat
  private plot_particles

  public Time_evolution

contains

  subroutine Time_evolution()

    integer :: i, flag
    
    flag = 1


    call relation_part(1)

    do i = 1, number_timesteps
       call Verlert_alg(i)
       vel_sqr(i) = sum(velocity**2)
       call thermostat(i)
       call plot_particles()
    end do
    
  end subroutine Time_evolution

  subroutine Verlert_alg(i)

    integer, intent(in) :: i

    velocity = velocity + forces*time_step/2.0
    position = modulo(position+velocity*time_step, length)
    call relation_part(i)
    velocity = velocity + forces*time_step/2.0
    
  end subroutine Verlert_alg

  subroutine thermostat(steps)

    integer, intent(in) :: steps

    if(steps < time_cut) then
       if(modulo(steps,20) == 0) then    
          velocity = velocity*sqrt((3*(num_particles-1)*temp_target/vel_sqr(steps)))
       end if
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
