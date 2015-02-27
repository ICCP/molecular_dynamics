program ArgonGas

  use global
  use initialization
  use time_evol
  use energy
  use pressure
  use pair_correl

  implicit none

  !$$ First we initilize the problem
  call Initialize_problem() 

  !$$ Second we start the simulation and let the system
  !$$ reach equilibrium
  call Time_evolution()

  !$$ Third we continue the simulation and determine
  !$$ physical quantities
  call write_files()
  call print_pair_corre()
  call calc_print_mean_press()
  call calc_print_heat_capac()

!!$  call EndPlot()

end program ArgonGas
