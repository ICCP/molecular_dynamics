module const_temp

  implicit none

  private

  public thermostat

contains

  subroutine thermostat(particles, temp_target, velocity, steps)
    integer, intent(in) :: particles
    integer, intent(in) :: steps
    real(8), intent(in) :: temp_target
    real(8), intent(inout) :: velocity(3,particles)


    real(8) :: temp_final

    temp_final = sum(velocity**2)/(3*(particles-1))

    if(modulo(steps,40) == 0 .and. steps < 700) then    
       velocity = velocity*sqrt(temp_target/temp_final)
    end if

    write (14,*) steps, temp_final

    
  end subroutine thermostat
  
end module const_temp
