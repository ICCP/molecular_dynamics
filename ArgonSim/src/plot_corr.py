#!/usr/bin/env python
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
import matplotlib.patches as mpatches

SaveFigure = True #True will save as png instead of displaying

#def plotaxis(large): 
 # x = np.linspace(-large, large, 100)
 # y = x*0.0
 # z = x*0.0
  #   ax.plot(x, y, z)
  #y = np.linspace(-large, large, 100)
 #x = x*0.0
 # z = x*0.0
  #   ax.plot(x, y, z)
 # z = np.linspace(-large, large, 100)
 # x = x*0.0
 # y = x*0.0
  #   ax.plot(x, y, z)

for desiredstart in [0,1,2,3]:
  fig, ax1 = plt.subplots()

  f = open("data/pair_correlation.dat")
  i = 1
  Data = []
  timeData = []
  ctr = 0
  current = 0
  for line in f:
    if (line[1] == "X"):
      print ("hello")
      current += 1
    else:
      if (current == desiredstart and ctr < 50):
        Data.append(line)
        timeData.append(10*ctr)
        ctr = ctr+1
    
  lns1 = ax1.plot(timeData, Data)
  x1,x2,y1,y2 = plt.axis()
  plt.axis((x1,x2,0,.05))
  ax1.set_xlabel('Distance (Sigmas)')
  ax1.set_ylabel('Number of Occurances')


  #plt.show()

  f.close()

  if SaveFigure:
    plt.savefig('plots/correlation'+str(desiredstart)+'.png', bbox_inches='tight')


  plt.close()

