program ArgonGas

  use global
  use initialization
  use time_evol
  use force
  use energy
  use pressure
  use pair_correl

  implicit none

  call Initialize_problem() 
    
  call Time_evolution()

!!$  call EndPlot()

end program ArgonGas
