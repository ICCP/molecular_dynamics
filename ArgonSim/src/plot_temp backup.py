#!/usr/bin/env python
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

SaveFigure = False #True will save as png instead of displaying
def plotaxis(large): 
  x = np.linspace(-large, large, 100)
  y = x*0.0
  z = x*0.0
  #   ax.plot(x, y, z)
  y = np.linspace(-large, large, 100)
  x = x*0.0
  z = x*0.0
  #   ax.plot(x, y, z)
  z = np.linspace(-large, large, 100)
  x = x*0.0
  y = x*0.0
  #   ax.plot(x, y, z)


f = open("data/temp.dat")
data = []
for line in f:
  data.append(line)

plt.plot(data)
plt.show()

f.close()

if SaveFigure:
  plt.savefig('data/temp.png', bbox_inches='tight')

plt.close()

