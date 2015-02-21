program ArgonGas

  use global
  use initialization
  use time_evol
  use energy
  use pressure
  use pair_correl

  implicit none
 
  call Initialize_problem() 
  call Time_evolution()

  call write_files()
  call print_pair_corre()
  call calc_print_mean_press()
  call calc_print_heat_capac()

!!$  call EndPlot()

end program ArgonGas
