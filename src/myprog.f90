program ArgonGas

  use initialization
  use force
  use energy
  use const_temp

  implicit none

  integer, parameter :: particles = 864
  real(8), parameter :: density = 0.8
  real(8), parameter :: temp_target = 1

  integer, parameter :: boxes = nint((particles/4)**(1.0/3))
  integer, parameter :: number_timesteps = 2000
  real(8), parameter :: time_step = 0.004
  real(8), parameter :: length = boxes*(4.0/density)**(1.0/3)
  real(8), parameter :: PI = 4*atan(1.0)

  real(8) :: position(3,particles), velocity(3,particles), forces(3,particles)
  real(8) :: pair_corre(nint(length/2*100))
  real(8) :: ener_kin, ener_pot

  integer, parameter :: time_sections = 325
  real(8) :: pressure(number_timesteps), mp, sdp
  real(8) :: mean_preassure(int((number_timesteps-700)/time_sections))
  real(8) :: std_dev(int((number_timesteps-700)/time_sections))

  integer :: i,j, flag

  print* , int((number_timesteps-700)/time_sections), (number_timesteps-700.0)/time_sections
  flag = 1
    
  call initialize_position(position, particles, density)
  call initialize_velocity(velocity, particles)
  call setting_cero_velocity(velocity, particles)
  forces = 0._8
  pair_corre = 0d0

  open (unit=11,file='ener_kin_data.txt')
  open (unit=12,file='ener_pot_data.txt')
  open (unit=13,file='ener_tot_data.txt')
  open (unit=14,file='temp_final_data.txt')
  open (unit=15,file='pair_corre_data.txt')
  open (unit=16,file='pressure_data.txt')

  call initplot ('lightblue', 800,800, 'out.ps', 1)
  call Framing (0._8, 0._8, length, length)
  call putstopbutton()

  call calculate_force(forces, particles, position, length, ener_pot, pair_corre, flag, pressure(1))

  do i = 1, number_timesteps
     velocity = velocity + forces*time_step/2.0
     position = modulo(position+velocity*time_step, length)
     call calculate_force(forces, particles, position, length, ener_pot, pair_corre, flag, pressure(i))
     velocity = velocity + forces*time_step/2.0
     call plot_particles(position, particles)
     ener_kin = sum(velocity**2*0.5)
     call thermostat(particles, temp_target, velocity, i)

     pressure(i) = pressure(i)/(3*particles*length**3)
     pressure(i) = density*temp_target - pressure(i)
     pressure(i) = pressure(i) + 16/3*PI*density**2*(2/3*3.2**(-9)-3.2**(-3))

     if (i == 700) then
        flag = 0
     end if

     write(11,*) i, ener_kin
     write(12,*) i, ener_pot
     write(13,*) i, ener_kin+ener_pot
     write(16,*) i, pressure(i)
  end do
 
  pair_corre = pair_corre*2.0*(length**3)
  pair_corre = pair_corre/(particles*(particles-1))
  pair_corre = pair_corre/(number_timesteps-700)

  do j = 1, nint(length*0.5*100)
     write(15,*) j/100.0, pair_corre(j)
  end do

  !calculo de la presion promedio y de la desviacion estandar
  mean_preassure = 0d0
  std_dev = 0d0
  i=1
  mp = 0.0
  sdp = 0.0
  do j = 700,  number_timesteps
     mean_preassure(i) = mean_preassure(i) + pressure(j)
     std_dev(i) = std_dev(i) + pressure(j)**2
     if (modulo(j-700+1, time_sections) == 0) then
        mean_preassure(i) = mean_preassure(i) / time_sections
        std_dev(i) = std_dev(i) / time_sections
        std_dev(i) = std_dev(i)-mean_preassure(i)**2
        mp = mp + mean_preassure(i)
        sdp = sdp + sqrt(std_dev(i))
        print*, j, mean_preassure(i), sqrt(std_dev(i))
        i = i+1
     end if
  end do

  mp = mp/int((number_timesteps-700)/time_sections)
  sdp = sdp/sqrt(1.0*int((number_timesteps-700)/time_sections))

  print*, mp, sqrt(sdp)
  
  call EndPlot()



contains
  
  subroutine plot_particles(pos, n)
    integer :: i, n
    real(8) :: pos(3,n)
    
    do i= 1, n
       call setpoint(pos(1,i), pos(2,i))
    end do
  end subroutine plot_particles

end program ArgonGas
