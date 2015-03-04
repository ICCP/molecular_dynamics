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

  f = open("data/energy.dat")
  i = 1
  KData = []
  PData = []
  TotData = []
  timeData = []
  ctr = 0
  current = 0
  for line in f:
    if (line[3] == "X"):
      print ("hello")
      current += 1
    else:
      if (current == desiredstart):
        if (ctr > 180):
          if ((ctr%3)==0):
            KData.append(line)
            timeData.append(8.65*int(ctr/3))
          if ((ctr%3)==1):
            PData.append(line)
          if ((ctr%3)==2):
            TotData.append(line)
        ctr = ctr+1
    
  lns1 = ax1.plot(timeData, KData, color = 'g', label = "Kinetic Energy")
  timeData.append(.004+timeData[-1])
  ax2 = ax1.twinx()
  lns2 = ax2.plot( timeData, PData, color = "r", label = "Potential Energy")
  lns3 = ax2.plot( timeData, TotData,color = 'b', label = "Total Energy")
  ax1.set_xlabel('Time (fs)')
  ax1.set_ylabel('Kinetic Energy (Natural Units)', color='g')
  ax2.set_ylabel('Potential & Total Energy (Natural Units)', color='r')

  lns = lns1+lns2+lns3
  labs = [l.get_label() for l in lns]
  ax1.legend(lns, labs, loc=4)

  #plt.show()

  f.close()

  if SaveFigure:
    plt.savefig('plots/energy'+str(desiredstart)+'.png', bbox_inches='tight')


  plt.close()

