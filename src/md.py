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
# Assign values from array
time  = energy_arr[:,0]
PE = energy_arr[:,1]
KE = energy_arr[:,2]
tot_E = energy_arr[:,3]
Temp = energy_arr[:,4]
momentum = energy_arr[:,5]
#
#plt.figure(1,figsize=(6,6))
plt.plot(time,PE,'r',label='Potential Energy')
plt.plot(time,KE,'b',label='Kinetic Energy')
plt.plot(time,tot_E,'g',label='Total Energy')
plt.legend(loc=0)
plt.xlabel('time')
plt.ylabel('Energy')
plt.title('Energy vs time plot')
plt.show()
#plt.savefig("Energyplot.pdf")
#
#plt.figure(1,figsize=(6,6))
plt.plot(time,Temp)
plt.xlabel('time')
plt.ylabel('Scaled Temperature')
plt.title('Temperature vs time plot')
plt.show()
#plt.savefig("TempPlot.pdf")
#
#plt.figure(1,figsize=(10,6))
plt.plot(time,momentum)
plt.xlabel('time')
plt.ylabel('Momentum of Center of Mass')
plt.title('Momentum vs time plot')
plt.show()
#plt.savefig("MomPlot.pdf")
#
with open('Velocity.dat', 'r') as vel_file:
    velocity_data =[]
    for line in vel_file:
        line = line.strip()
        columns = line.split()
        velocity_data.append(columns)
vel_data = np.array(velocity_data)
x  = vel_data[:,0]
v1 = vel_data[:,1]
v2 = vel_data[:,2]
v3 = vel_data[:,3]
v4 = vel_data[:,4]
v5 = vel_data[:,5]
#plt.figure(1,figsize=(6,6))
plt.plot(x,v1,'r',x,v2,'b',x,v3,'k',x,v4,'g',x,v5,'m')
plt.xlabel('time')
plt.ylabel('Velocity of individual particle')
plt.title('Velocity vs time plot')
plt.show()
#plt.savefig("VelPlot.pdf")
