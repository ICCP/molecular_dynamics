program ArgonGas

  use global
  use initialization
  use time_evol
  use energy
  use pressure
  use pair_correl

  implicit none

  real(8) :: mean_pressure, error_pressure
  
  call Initialize_problem() 
  call Time_evolution()

  call write_files()
  call print_pair_corre()
  call mean_and_std_dev(pressures, 300, mean_pressure, error_pressure)
  call calc_print_mean_press(mean_pressure, error_pressure)

  call EndPlot()

end program ArgonGas
