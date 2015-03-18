#!/usr/bin/python
import sys
# from pylab import *
import numpy as np
import matplotlib.pyplot as plt
# This program performs MD simulation of N atoms in a 3-Dimensional space
# Program written by Tridip Das (dastridi)
# This program calculates magnetization with temperature variation 
# Susceptibility, Energy and Specific heat with temperature
# Also finds the critical temperature at which magnetization is
# lost. Scaled temperature is used, Boltzmann const, kB assumed as 1
m = 1
kB = 1
sigma = 1
epsilon  = 1
# Lattice position
init_lat_pos = open('initial_lattice.dat', 'w')
# Provide a number of latticein x direction below:
L = 2 
B = 2**(1/6.0)*sigma
Lattice = (2**0.5)*B # Length of a single lattice
#
print "Lattice ", Lattice
N = 4*L**3 # Total number of atoms
#
atom_pos = np.zeros((N,3))
# 

C = -1 # Counter for every atom
for i in range(0,L):
	for j in range(0,L):
		for k in range(0,L):
			arr1 = np.array([i, j, k]) # To store lattice corner			

			C = C+1  
			atom_pos[C,:] = Lattice * arr1
			atom_pos[C+1,:] = Lattice * (arr1 + np.array([0.5, 0.5, 0.0]))
			atom_pos[C+2,:] = Lattice * (arr1 + np.array([0.5, 0.0, 0.5]))
			atom_pos[C+3,:] = Lattice * (arr1 + np.array([0.0, 0.5, 0.5]))
			C = C+3

# print """ \n""", atom_pos
np.savetxt(init_lat_pos,atom_pos,fmt='%f')
init_lat_pos.close()
#
try:  
   with open('Energy.dat', 'r') as datafile:
	first_line = datafile.readline()
	energy_data =[]
	for line in datafile:
		line = line.strip()
		columns = line.split()	 
		energy_data.append(columns)
   datafile.close()
except:
   print " Please RUN mdfort.exe to generate Energy.dat file "
   print " Once Energy.dat file created rerun md.py "
   exit()
#
energy_arr = np.array(energy_data)
#
x  = energy_arr[:,0]
y1 = energy_arr[:,1]
y2 = energy_arr[:,2]
y3 = energy_arr[:,3]
y4 = energy_arr[:,4]
y5 = energy_arr[:,5]
plt.plot(x,y1,'r',label='PE')
plt.plot(x,y2,'b',label='KE')
plt.plot(x,y3,'g',label='E')
plt.plot(x,y4,'c',label='T')
# plt.plot(x,y1,'r',x,y2,'b',x,y3,'g',x,y4,'c')
plt.legend(loc=0)
plt.xlabel('time')
plt.ylabel('Energy')
plt.title('Energy vs time plot');
plt.savefig("Energyplot.pdf")
plt.show()
