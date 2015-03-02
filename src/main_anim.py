#IMPORTS
import math
import numpy as np
import f90force_corr
import initial_conditions as ic
import plot
import get
import anim_md

#PARAMETERS
Ttarget=0.5
M = 1
N=4*M ** 3 # Choose from 4, 32, 108, 256, 500, 864
dt=0.004
simulation_time=1
m=1
sigma=1
eps=1
rc=10000
r0=sigma*math.pow(2.0,1.0/6.0)
kb=1
n=2*M
l = r0 / (math.pow(2.0, 1.0/2.0)) # distance between two "lattice points"
L=l*n
V = L ** 3.0
steps = 1000
numSteps = int(simulation_time/dt)

correlation = np.zeros((steps),dtype=float)
forces = np.zeros((N,3),dtype=float)
potential = 0.0
totE = np.zeros(numSteps) 
KE = np.zeros(numSteps)
PE = np.zeros(numSteps)
x = np.zeros(numSteps)
T = np.zeros(numSteps)

#ALGORYTHM

#Initial Conditions
positions = ic.set_initial_positions(N, sigma, l, n, L)
momenta = ic.set_initial_momenta(N, kb, Ttarget, m)
forces,correlation,potential = f90force_corr.update_forces(positions,forces,correlation,potential,rc,L,[N,steps])


     
def simulate(N, L, positions, momenta, dt, potential, rc, steps, forces, correlation, Ttarget):
    
	#integrator
	forces,correlation,potential = f90force_corr.update_forces(positions,forces,correlation,potential,rc,L,[N,steps])
	momenta += 0.5 * dt * forces
	positions = (positions + dt*momenta/m ) % L
	forces,correlation,potential = f90force_corr.update_forces(positions,forces,correlation,potential,rc,L,[N,steps])
	momenta += 0.5 * dt * forces

	#thermostat
	T = m* np.sum(momenta ** 2) / ((N-1) * 3 * kb)
	momenta *= math.sqrt(Ttarget / T)

	return positions, momenta

anim = anim_md.AnimatedScatter(N, L, positions, momenta, simulate, dt, potential, rc, steps, forces, correlation, Ttarget)
anim.show()


 
     
    
    




