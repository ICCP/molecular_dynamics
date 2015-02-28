#Initial momenta module
import numpy as np
import random as rn

def init_mom(N, temp):
  momenta = np.random.normal(0.0, np.sqrt(temp), (N,3))
  avmom = np.zeros((1,3),dtype=float)
  avmom = np.sum(momenta, axis=0)/N
  momenta = momenta - avmom
  return momenta

def init_pos(N,density):
  npdim = int(round((N/4)**(1./3.)))
#  print npdim, 'npdim'
  l = (N/density)**(1.0/3)
  a = 0.5*l/npdim
#  print npdim, l, a
  pos = np.zeros((N,3), dtype = float)
  cnt = 0
  for k in range(npdim):
    for j in range(npdim):
      for i in range(npdim):
        pos[cnt,:]   = [0+2*a*i,0+2*a*j,0+2*a*k]
        pos[cnt+1,:] = [a+2*a*i,a+2*a*j,0+2*a*k]
        pos[cnt+2,:] = [a+2*a*i,0+2*a*j,a+2*a*k]
        pos[cnt+3,:] = [0+2*a*i,a+2*a*j,a+2*a*k]
        cnt = cnt + 4

  return pos, l, npdim

def init_forc(N):
  forces = np.zeros((N,3), dtype = float)
  return forces

def init_dist(N):
  distances = np.zeros((N,N), dtype = float)
  return distances

def init_bins(nbins,lpnum):
  bin_vec_tot = np.zeros(nbins, dtype = float)
  finalbins = np.zeros(nbins, dtype = float)
  return bin_vec_tot,finalbins
  
def init_toten(lpnum):
  toten = np.zeros((lpnum,1), dtype = float)
  return toten
  
def init_presvirialtime(lpnum):
  presvirial = np.zeros((lpnum,1), dtype = float)
  return presvirial
  
def init_kenarray(lpnum):
  kenarray = np.zeros((lpnum,1), dtype = float)
  return kenarray
  
def init_potarray(lpnum):
  potarray = np.zeros((lpnum,1), dtype = float)
  return potarray
