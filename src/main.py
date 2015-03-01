####### IMPORTS
import math
import numpy as np
import f90force as f90
# import forces.f90 as f90
# import forces2 as f90
import initial_conditions as ic
import plot
import get

####### PARAMETERS
# Constants
Ttarget = 5				# Temperature maintained with thermostat
M = 6					# Number of unit cells in the length of the box
N = 4*M ** 3  				# Number of particles in the box
dt = 0.001				# Time step, ideal value is 0.004
simulation_time = 1			# Simulation time
m = 1					# Mass of the particles
sigma = 1				# Sigma parameter for the Lennard-Jones potential
eps = 1					# Epsilon parameter for the Lennard-Jones potential
r0 = sigma * math.pow(2.0,1.0/6.0)	# Initial distance between the particles
kb = 1					# Boltzmann constant
n = 2*M					# Number of "lattice points"
l = r0 / (math.pow(2.0, 1.0/2.0)) 	# Distance between two "lattice points"
L = l * n				# Length of one side of the box
rc= 2.5					# Distance under which we start to calculate force and study pair-correlation
V = L ** 3.0				# Volume of the box
correlation_steps = 1000		# Number of shells used in the space discretization when calculating correlation
numSteps = int(simulation_time/dt)	# Number of time steps

# Variables
correlation = np.zeros(correlation_steps)	# Correlation function
forces = np.zeros((N,3))			# Vector containing the force vectors for each particle
PE = np.zeros(numSteps)				# Total potential in the box
totE = np.zeros(numSteps) 			# Total energy
KE = np.zeros(numSteps)				# Total kinetic energy (sum over all partcles of individual kinetic energies)
x = np.zeros(numSteps)				# Vector of times to plot energies
T = np.zeros(numSteps)				# Temperature the box
P = np.zeros(numSteps)				# Pressure the box

####### ALGORYTHM

#Initial Conditions
positions = ic.set_initial_positions(N, r0, l, n, L)	# The vector of individual positions is initialized in a fcc grid
momenta = ic.set_initial_momenta(N, kb, Ttarget, m)	# The vector of individual momenta is initialized according to a Maxwell distribution
forces,correlation,PE[0],P[0] = f90.update_forces(positions,forces,correlation,0,rc,L,0,kb,T,V,[N,correlation_steps])
							# Force, correlation and potential are initialized with the initial position/momenta

     
for i in xrange(numSteps):
    # Calculation of relevant quantities
    T[i] = m* np.sum(momenta ** 2) / ((N-1) * 3 * kb)
    KE[i] = np.sum(momenta ** 2 ) /(2.0*m)
    x[i] = i*dt


    # Integrator (Verlet algorithm)
    momenta += 0.5 * dt * forces
    positions = (positions + dt*momenta/m ) % L
    forces,correlation,PE[i],P[i] = f90.update_forces(positions,forces,correlation,0,rc,L,0,kb,T[i],V,[N,correlation_steps])
    momenta += 0.5 * dt * forces
   
   
    # Thermostat maintains T = Ttarget initially
    if ((i % 10 == 0)&(i*dt < 0.4)) :
        momenta *= math.sqrt(Ttarget / T[i])

# Results
print 'Mean value for temperature: ',get.mean(T),' with variance of: ', get.var(T)
print 'Mean value for pressure: ',get.mean(P),' with variance of: ', get.var(P)
plot.correlation_function(correlation, correlation_steps, N, L, numSteps)
plot.vector(x,KE+PE)
plot.vector(x,T)

