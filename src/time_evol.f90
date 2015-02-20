module time_evol

  use global
  use force
  use energy
  use pressure
  use pair_correl

  implicit none

  private thermostat
!!$  private plot_particles

  public Time_evolution

contains

  subroutine Time_evolution()

    real(8) :: pair_corre(nint(length/2*100))
    real(8) :: pressure(number_timesteps)
    
    real(8) :: ener_pot
    real(8) :: mean_pressure, stdev_pressure
    integer :: i, flag
    
    flag = 1
    pair_corre = 0d0
    call relation_part(ener_pot, pair_corre, flag, pressure(1))

    do i = 1, number_timesteps
       velocity = velocity + forces*time_step/2.0
       position = modulo(position+velocity*time_step, length)
       call relation_part(ener_pot, pair_corre, flag, pressure(i))
       velocity = velocity + forces*time_step/2.0
       
       call write_energy(ener_pot, i)
       call thermostat(i, flag)
       call print_press(pressure(i), i)
!!$       call plot_particles()
    end do
    
    call print_pair_corre(pair_corre)
    
    call mean_and_std_dev(pressure, 300, mean_pressure, stdev_pressure)

    mean_pressure = mean_pressure/(3*num_particles*temp_target)
    mean_pressure = 1 + mean_pressure
    mean_pressure = mean_pressure + 16*PI*num_particles/(3*length**3*temp_target)*(2/3*3.2**(-9)-3.2**(-3))
    stdev_pressure = stdev_pressure/(3*num_particles*1.095)

    print*, "Average preassure:", mean_pressure
    print*, "Standar error of mean:", stdev_pressure
    
  end subroutine Time_evolution

  subroutine thermostat(steps, flag)

    integer, intent(in) :: steps
    integer, intent(out) :: flag


    real(8) :: temp_final

    temp_final = sum(velocity**2)/(3*(num_particles-1))
    flag = 0

    if(steps < 700) then
       flag = 1
       if(modulo(steps,40) == 0) then    
          velocity = velocity*sqrt(temp_target/temp_final)
       end if
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
