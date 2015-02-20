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

    
    real(8) :: mean_pressure, stdev_pressure
    integer :: i, flag
    
    flag = 1
    call relation_part(1)

    do i = 1, number_timesteps
       call Verlert_alg(i)
       kin_energy(i) = sum(velocity**2*0.5)
       call thermostat(i)
!!$       call plot_particles()
    end do
    
    call write_energy_pressure()
 
    call print_pair_corre()
    
    call mean_and_std_dev(pressures, 300, mean_pressure, stdev_pressure)

    mean_pressure = mean_pressure/(3*num_particles*temp_target)
    mean_pressure = 1 + mean_pressure
    mean_pressure = mean_pressure + 16*PI*num_particles/(3*length**3*temp_target)*(2/3*3.2**(-9)-3.2**(-3))
    stdev_pressure = stdev_pressure/(3*num_particles*1.095)

    print*, "Average preassure:", mean_pressure
    print*, "Standar error of mean:", stdev_pressure
    
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

    real(8) :: temp_final

    temp_final = sum(velocity**2)/(3*(num_particles-1))
 
    if(steps < 700) then
       if(modulo(steps,40) == 0) then    
          velocity = velocity*sqrt(temp_target/temp_final)
       end if
    else if(steps == 700) then
       flag = 0
    end if


    write (14,*) steps, temp_final
    
  end subroutine thermostat

!!$  subroutine plot_particles()
!!$
!!$    integer :: i
!!$
!!$    do i= 1, num_particles
!!$       call setpoint(position(1,i), position(2,i))
!!$    end do
!!$
!!$  end subroutine plot_particles
  
end module time_evol
