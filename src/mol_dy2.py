# Main file

# import python classes

import numpy as np
import random as rn
import math
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import sys

# import self produced classes

import forcemodule as fm
import init_sys
import correlationmodule
import pressuremodule

# flags

corrflag = 2 # 2 to obtain pressure, 1 to obtain correlation function

# independent parameters

dt = 0.004
N=864
lpnum = 10000
cutoff = int(lpnum/5.)
density = 1.2
temp = 0.8
Ttarg = 0.5
Kb =  1
nbins = 250 #number of radial shells
n = 100 #number of timesteps per timeblock in pressure calculation
presint = 0.97596 #This is the value for the integral which is the third term to obtain the pressure


# Loading initial conditions

mom = init_sys.init_mom(N, temp) 
pos, l, npdim = init_sys.init_pos(N, density) 
forces = init_sys.init_forc(N)
distances = init_sys.init_dist(N)
bin_vec_tot, finalbins = init_sys.init_bins(nbins,lpnum)
toten = init_sys.init_toten(lpnum)
prestime = init_sys.init_presvirialtime(lpnum-cutoff)
kenarray = init_sys.init_kenarray(lpnum)
potarray = init_sys.init_potarray(lpnum)

# Iteration Verlet method

# Initializing the first values for the loop
pot = 0.0
presvirial = 0.0
forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances,presvirial,[N])
lamda=0
cnt = 0

for lp in range(lpnum):
  mom = mom + forces*0.5*dt
  pos = (pos + mom*dt) % l           # % l means modulo of l, hence it adds/subtracts n*l untill 0<pos<l
  forces, pot, distances, presvirial = fm.calc_forces(pos,forces,pot,l,distances, presvirial, [N])
  mom = mom + forces*0.5*dt
  Ken = np.sum(mom*mom*0.5)
  kenarray[lp,:] = Ken
  potarray[lp,:] = pot
  toten[lp,:] = Ken + pot
  if lp < (cutoff):
    if lp % 10 == 0:
      lamda = np.sqrt((Ttarg*3.*(N-1.)*Kb)/(np.sum(Ken)*2.))
      mom = mom*lamda
  elif corrflag == 1:
    corrflag = 0
    bin_vec_tot = bin_vec_tot + correlationmodule.cor(npdim,N,distances,nbins,finalbins,corrflag,Ttarg,density)
    corrflag = 1
  elif corrflag == 2:
    prestime[cnt,:] = density*((2*Ken/(3*N)) + presvirial/(3*N) + np.pi*N*presint/Ken) 
    cnt = cnt + 1

# pressure module, writes data into a datafile

if corrflag == 2:
  pressuremodule.pres(lpnum,cutoff,n,prestime,toten,kenarray,potarray,Ttarg,density,N)

# Calculate and plot correlationfunction

if corrflag == 1:
  finalbins = 2*bin_vec_tot/((lpnum - cutoff)*density*(N-1))
  bin_vec = correlationmodule.cor(npdim,N,distances,nbins,finalbins,corrflag,Ttarg,density)





















