#!/usr/bin/env python
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

SaveFigure = True #True will save as png instead of displaying
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


f = open("data/pair_correlation.dat")
Data = []
DataBackUp = []
i = 0
for line in f:
    Data.append(int(line))
    if (int(line) > 0):
      DataBackUp = Data
      latest = i
    i += 1

f.close()


# the histogram of the data
plt.plot(DataBackUp)
plt.axis([0,latest*1.25,0,max(DataBackUp)*1.25])
#plt.xlabel('dr')
#plt.ylabel('frequency')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([0,max(Data),0,length(Data)])
#plt.grid(True)

plt.show()


if SaveFigure:
  plt.savefig('data/corr.png', bbox_inches='tight')


plt.close()

