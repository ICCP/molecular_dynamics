import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

# Plots inputed vector x
def vector(x,y):
    plt.plot(x,y,hold='false')
    plt.show()

# Plots the position of the particles in 2D
def position(positions, L):  
    fig = plt.figure()
    ax = fig.gca() 
    for i in xrange(len(positions)):
        x, y, z = positions[i]
        ax.scatter(x,y)
    plt.xlim([0,L])
    plt.ylim([0,L])
    plt.show()

# Applies a normalization factor to each component of the correlation vector, then plots the correlation
def correlation_function(correlation, steps, N, L, numSteps):    
    
    max_distance = L / 2	# Largest distance between particles condidered in the calculation of the correlation function
    dr = max_distance / steps	# Width of each shell
    x_axis = np.zeros(steps)	# x-axis in the Correlation function plot
    for k in xrange(1,steps):
        r = k * dr
	correlation[k] *= 2*L**3/(N*(N-1)*4*math.pi*r**2*dr*(numSteps))
        x_axis[k] = r
     
    plt.bar(x_axis,correlation, width = max_distance / steps, bottom = 0)
    plt.show()
