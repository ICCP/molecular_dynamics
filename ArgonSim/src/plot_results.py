#!/usr/bin/env python
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

SaveFigure = True #True will save as png instead of displaying
numFigs = 100 #More than one will input/output multiple datafiles (use only with SaveFigure==True)

for n in range(0,numFigs):

  def plotaxis(large): 
     x = np.linspace(-large, large, 100)
     y = x*0.0
     z = x*0.0
     y = np.linspace(-large, large, 100)
     x = x*0.0
     z = x*0.0
     z = np.linspace(-large, large, 100)
     x = x*0.0
     y = x*0.0

  currOrder = 0
  while ((n+1) > ((10**currOrder)-1)):
    currOrder = currOrder + 1
  numSpaces = 5-currOrder

  f = open("data/position"+(" "*numSpaces)+str(n+1)+".dat")
  i = 1
  data = []
  for line in f:
    data_aux = []
    for x in line.split():
        data_aux.append(float(x))
    data.append(data_aux)
  f.close()
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  plotaxis(30)
  P.xlim([0,3.2])
  P.ylim([0,3.2])
  ax.set_zlim(0,3.2)
  if (n%5==0):
    print("Plot image creation: step ",n," out of ",numFigs)
  for point in data:
    ax.scatter(point[0], point[1], point[2])

  if SaveFigure:
    plt.savefig('data/PositionData'+("0"*numSpaces)+str(n+1)+'.png', bbox_inches='tight')
  else:
    plt.show()

  plt.close()

