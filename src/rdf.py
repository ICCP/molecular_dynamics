#!/usr/bin/python
import sys
# from pylab import *
import numpy as np
from math import sqrt, ceil, floor, pi
import matplotlib.pyplot as plt
# This program calculates pair distribution function (PDF) of N atoms in a 3-Dimensional space with periodic boundary condition
# Program written by Tridip Das (dastridi)
def calculateRDF(datafile):
 sigma = 1.0
 B = 2**(1/6.0)*sigma
 Lattice = (2**0.5)*B # Length of a single lattice
# Provide a number of lattice in x direction below:
 L = 6
 N = 4*L**3 # Total number of atoms
 Volume = (L*Lattice)**3
 density = N/Volume
#
 with open(datafile, 'r') as lat_file:
    lattice_data =[]
    for line in lat_file:
        line = line.strip()
        columns = line.split()
        lattice_data.append(columns)
 lat_file.close()
 atom_pos = np.array(lattice_data, dtype="float")
#
 side = L * Lattice
 half_side = side/2.0
 dr = 0.1
 max_bin = int(ceil(sqrt(3)*half_side/dr))
 print max_bin
#
 hist = np.zeros(max_bin+1)
 for i in range (N-1):
	for j in range(i+1, N):
		Rx = atom_pos[i,0] - atom_pos[j,0]
		Ry = atom_pos[i,1] - atom_pos[j,1]
		Rz = atom_pos[i,2] - atom_pos[j,2]		
		#
		# minimum image convention
		if ( Rx < - half_side): Rx = Rx + half_side
		if ( Rx  >  half_side): Rx = Rx - half_side
		if ( Ry < - half_side): Ry = Ry + half_side
                if ( Ry  >  half_side): Ry = Ry - half_side
		if ( Rz < - half_side): Rz = Rz + half_side
                if ( Rz  >  half_side): Rz = Rz - half_side
		# distance between particles
		Rij = sqrt(Rx*Rx + Ry*Ry + Rz*Rz)
		bin = int(ceil(Rij/dr)) # determine in which bin particle is
		hist[bin] += 1
 pos = np.zeros(max_bin+1)
# Normalize histogram data
 for i in range(1, max_bin+1):
	vol = ((i+1)**3 - i**3)*dr**3
	shell_vol = (4/3.0)*3.1416*vol
	hist[i] = hist[i]/(N*shell_vol*density)
	pos[i] = i*dr
 plt.plot(pos,hist)
 plt.xlabel('radial position')
 plt.ylabel('g(r)')
 plt.title('Radial Distribution Function plot')
 plt.show()
# plt.savefig("RDFPlot.pdf")
# run for different files
rdf = calculateRDF('final_lat.dat')
