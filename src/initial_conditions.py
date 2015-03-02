import math
import numpy as np
import random

# We initialize the postions of the particles in a fcc grid and return a vector of the positions of every individual particle. The particles are numbered from 0 to N-1, form left to right, top to bottom, and then front to back
def set_initial_positions(N, r0, l,n, L):
    
    positions = []     
    for k in xrange(n):		# We go throught the n "lattice points" in the z direction, labelled by k
        for j in xrange(n):	# We go throught the n "lattice points" in the y direction, labelled by j
            if (j+k) % 2 == 0:	# If j+k is even, the atoms are placed every two lattice points, starting at the first point
                for i in xrange(n/2):
                    positions.append([i*2*l,j*l,k*l])
            else:		# If j+k is odd, the atoms are placed every two lattice points, starting at the second point
                for i in xrange(n/2):
                    positions.append([i*2*l+l,j*l,k*l])
                    
    positions = np.array(positions) 
    return positions



# Returns a momenta for one particle according to a Maxwell distribution. The statistics vary with Boltzmanns constant, the temperature and the mass of the particle
def select_gaussian_momentum(kb, T, m):

    rand1=random.random() # Determines the sign of the momenta
    rand2=random.random() # Determines the absolute value of the momenta
    velocity = math.sqrt(-2 * kb * T * math.log(rand2) / m)
    if (rand1 > 0.5):
        return velocity
    else:
        return -velocity
    

# Returns a vector containing the individual momenta vectors for each particle chosen according to a Maxwell distribution. This function also corrects total momenta that could arise from the fact that the random choice of momenta is not evenly distributed.
def set_initial_momenta(N, kb, T, m):   
    momenta = []

# Construction of the vector
    for i in xrange(N):
        momentum = []
        for j in xrange(3):
            momentum.append(select_gaussian_momentum(kb, T, m))
            
        momenta.append(momentum)  
    
# Correction for initial momenta  
    for i in xrange(N):
        for j in xrange(3):
            momenta[i][j] -= np.sum(momenta,axis=0)[j]
    momenta = np.array(momenta)      
    return momenta
