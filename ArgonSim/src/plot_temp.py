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


fig, ax1 = plt.subplots()

f = open("data/temp.dat")
i = 1
Data = [[],[],[],[]]
timeData = []
ctr = 0
current = 0
for line in f:
  if (line[3] == "X"):
    current += 1
    ctr = 0
    if (current < 4):
      timeData = []
  else:
    if (ctr > 100):
      Data[current].append(line)
      timeData.append(8.65*ctr)
    ctr = ctr+1
for j in [1,2]:
  Data[j] = Data[j][:-9]
Data[0] = [0]*(len(timeData))
for i in [0,1,2,3]:
  print("Length of ",i,": ",len(Data[i]))
print("Length of time: ",len(timeData))  
lns1 = ax1.plot(timeData, Data[0], color = 'g', label = "T=0")
lns2 = ax1.plot( timeData, Data[1], color = "r", label = "T=1")
lns3 = ax1.plot( timeData, Data[2],color = 'b', label = "T=5")
lns4 = ax1.plot( timeData, Data[3],color = 'y', label = "T=10")
ax1.set_xlabel('Time (fs)')
ax1.set_ylabel('Temperature (Natural Units)')

lns = lns4+lns3+lns2+lns1
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=4)

#plt.show()

f.close()

if SaveFigure:
  plt.savefig('plots/temp.png', bbox_inches='tight')


plt.close()

