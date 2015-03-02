import numpy as np
import random
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.cm as cm
from scipy.stats import norm
import time

#Time the code
start_time = time.clock()

# Input parameters
global A, N, L, V, dens, T, m, Sig, Eps, t, state, stat, temp, press, Rc_inv

#Physical constants
dens = 1.0   #Density
T=0.2         #Temperature
m=1.0         #Particle mass
Sig=1        #Sigma in Lennard-Jones potential 0.3345e-9
Eps=1        #Epsilon in Lennard-Jones potential eps/kb 125.7 K

#System sizes
A = 7             #Number of FCC boxes
N = 4*A**3          #Number of atoms
L=(N/dens)**(1.0/3) #Length of the box in all three directions
V=L**3              #Box volume

#Options for the output
state = 0      #State of the system: 0 for animation, 1 for correlation
stat = 3       #State of plots:
               #0 for correlation,
               #1 for Temperature
               #2 for the pressure
               #3 for the Energy

#Computational parameters
t_final = 200   #Amount of steps of delta_t
Rc= L/5       #Cutoff radius #L/2.1   #should be about 5*Sig
Rc_inv=(Rc)**(-1) #Inverse of the cut-off raduis
delta_t=0.004  #Time-step
t_ren=20      #Number of steps before T-renormalization
maxi = t_final+1     #Maximum number of runs for plot of T and E

delta_r=0.02                        #Shell width for correlation function
b=(L/2.)/delta_r-((L/2.)/delta_r)%1 #Number of Correlation function bins
n_plot = 4                         # Number of correlation plots

temp = np.zeros((maxi,2),dtype=float)
press = np.zeros((maxi,2),dtype=float)
E_total = np.zeros((maxi,2),dtype=float)
E_kin = np.zeros((maxi,2),dtype=float)
E_pot = np.zeros((maxi,2),dtype=float)
############################################################

#Generate initial positions in FCC lattice for all N particles
def init_position(A, L):
  position = np.zeros((N,3),dtype=float)
  k = 0
  l = 0
  m = 0
  count = 0
  for k in range(A):
    for l in range(A):
      for m in range(A):
        position[count,:] = [k, l, m]
        position[count+1,:] = [k+0.5, l+0.5, m]
        position[count+2,:] = [k+0.5, l, m+0.5]
        position[count+3,:] = [k, l+0.5, m+0.5]
        count += 4
  position = position*L/A
  return position

############################################################

#Generate initial momenta with a gaussian distribution for all N particles
def init_momenta(T,N):
  mu=0
  sigma=(2*T)**(1.0/2)

  momenta = np.random.normal(mu,sigma,(N,3))

  avg = np.zeros((3),dtype=float)
  for i in range(0, 3):
    avg[i] = sum(momenta[:,i]/N)
    momenta[:,i] = [j - avg[i] for j in momenta[:,i]]
  return momenta

############################################################
########################CHECKS##############################
############################################################

# Momenta Conservation Check
def check_mom(momenta):
    total_mom=np.zeros((3),dtype=float)
    for x in range(3):
        total_mom[x]=sum(momenta[:,x])

    return total_mom

############################################################
####################ITERATIVE PROCESS#######################
############################################################

#Calculates the Lennard Jones Force for all particles
def calc_force(position, momenta,N,Sig,Eps):
  #Initializing all arrays and numbers, calling neccesary functions
  force_v=np.zeros((N,3),dtype=float)
  V=0.0
  p_vir=0
 
  mom_abssq=np.sum(momenta*momenta,axis=1)

  for i in range(1,N):
    dist=position[i,:]-position[:i,:]
    dist[dist>(L/2.0)]-=L
    dist[dist<-(L/2.0)]+=L

    dist_abs=np.sum(dist*dist,axis=1)
    dist_abs_inv=(dist_abs)**(-1)

    #implementing Cut-off radius
    dist_abs_inv[dist_abs_inv < Rc_inv]=0

    force = -24.0*Eps*(Sig**6.0*dist_abs_inv**(4.0)-2.0*Sig**12.0*dist_abs_inv**(7.0))
  
    force_temp=dist*force[:,np.newaxis]
    p_vir += np.sum(force_temp*dist)
    
    force_v[:i,:]-=force_temp
    force_v[i,:] = np.sum(force_temp,axis=0)
    
    V += np.sum(4.0*Eps*((Sig/dist_abs)**6.0-(Sig/dist_abs)**3.0))
  
  
  K=np.sum((0.5/m*mom_abssq),axis=0)
  total_E = V + K
  T_calc=(2.*K)/(3.*N)
  
  # Pressure Calculation with Virial Theorem
  pressure=L**(-3)*(N*T_calc+p_vir*(3.**(-1)))
  #print(pressure)  
  

  return force_v, V, total_E, K, pressure

############################################################
#autocorrelation function
def autocorr(x):
    mean = sum(x)/len(x)
    result = np.correlate(x-mean, x-mean, mode='full')
    return result[result.size/2:]

def standdev(x,tau,t_final):
    mean=sum(x)/len(x)
    mean2=sum(x*x)/len(x)
    stddev=(mean2-mean**2.)**0.5
    print(stddev)
    M=t_final*1.
    err=(2.*tau/M)**0.5*stddev
    return err, mean
####################################################################

#The correlation function
def cor_func(position,L,N):
  #Initialize arrays
  dist=np.zeros((3),dtype=float)
  dist_abs=np.zeros((N),dtype=float)
  #Choose one particle
  for i in range(N):
    dist=position[i]-position[0]
    #If particles are L/2 or more away, put them at the other side of the box.
    dist[dist>=(L/2.)]-=L
    dist[dist<-(L/2.)]+=L
    dist_abs[i]=sum(dist*dist)**1/2
  return dist_abs

#####################################################################

#Putting together all iterative steps
def iteration(position,momenta,force,initial_E, initial_V,delta_t,N,L,cnt):
    #print(cnt)

    #Verlet's Theorem
    momenta += 0.5*force*delta_t
    position += momenta*delta_t
      # Boundary conditions
    position=position%L
    force, V, total_E, K, pressure = calc_force(position,momenta,N,Sig,Eps)
    momenta += 0.5*force*delta_t
    
    #cnt += 1.
    # Saving the values for the plots at the end
    temp[cnt,1]=(2.*K)/(3.*N)
    temp[cnt,0]=cnt


    # Performing the renormalisation every t_ren time steps
    if cnt%t_ren==0:
        labda = (T/(temp[cnt,1]))**0.5
        momenta = momenta*labda
        #print(labda)

    #When the animation is not runned,
    #the program can calculate several physical properties
        
    if state==1:

      if cnt==t_final:
        print("Time(s):", time.clock() - start_time)
      if stat==0:   # Correlation function
          plot_nr=np.floor(t_final/n_plot)
          k=cnt%(plot_nr)        
          if k==0:
            hist=np.histogram(cor_func(position,L,N)[1:N],bins=b,density=True)
            xhist=hist[1]
            yhist=np.concatenate(([0],hist[0]),axis=1)/(4*np.pi*np.multiply(hist[1],hist[1]))
            fig=plt.plot(xhist,yhist)
            if cnt==t_final:
                plt.xlabel('Distance between the particles')
                plt.ylabel('Density')
                plt.title('Correlation Function')
                plt.show()
      elif stat==1: # Temperature
            if cnt==t_final:
                fig=plt.plot(temp[:,0],temp[:,1])
                plt.xlabel('Time [# of iteration steps]')
                plt.ylabel('Temperature [a.u.]')
                plt.title('Renormalized Temperature')
                plt.show()
                A=autocorr(temp[:,1])**2
                tau = np.where(A<=0.5*max(A))[0][0]
                err = standdev(temp[:,1],tau,t_final)[0]
                mean = standdev(temp[:,1],tau,t_final)[1]
                print("T is",mean,"error is",err)
      elif stat==2: # Pressure
            press[cnt,1]=pressure
            press[cnt,0]=cnt
            if cnt==t_final:
              fig=plt.plot(press[:,0],press[:,1])
              plt.xlabel('Time [# of iteration steps]')
              plt.ylabel('Pressure [a.u.]')
              plt.title('Pressure')
              plt.show()
              A=autocorr(press[:,1])**2
              tau = np.where(A<=0.5*max(A))[0][0]
              err = standdev(press[:,1],tau,t_final)[0]
              mean = standdev(press[:,1],tau,t_final)[1]
              print("P is",mean,"error is",err)
      elif stat==3: # Energies
            E_total[cnt,1]=total_E
            E_total[cnt,0]=cnt
            E_kin[cnt,1]=K
            E_kin[cnt,0]=cnt
            E_pot[cnt,1]=V
            E_pot[cnt,0]=cnt
            if cnt==t_final:
                tau=np.zeros((3,1),dtype=float)
                err=np.zeros((3,1),dtype=float)
                mean=np.zeros((3,1),dtype=float)
                fig=plt.plot(E_total[:,0],E_total[:,1],label='Total Energy')
                fig=plt.plot(E_kin[:,0],E_kin[:,1],label='Kinetic Energy')
                fig=plt.plot(E_pot[:,0],E_pot[:,1],label='Potential Energy')
                plt.xlabel('Time [# of iteration steps]')
                plt.ylabel('Energy [a.u.]')
                plt.title('Energy')
                plt.legend()
                plt.show()
                A=autocorr(E_total[:,1])**2
                tau[0] = np.where(A<=0.5*max(A))[0][0]
                err[0] = standdev(E_total[:,1],tau[0],t_final)[0]
                mean[0] = standdev(E_total[:,1],tau[0],t_final)[0]
                B=autocorr(E_kin[:,1])**2
                tau[1] = np.where(B<=0.5*max(B))[0][0]
                err[1] = standdev(E_kin[:,1],tau[1],t_final)[0]
                mean[1] = standdev(E_kin[:,1],tau[1],t_final)[1]
                C=autocorr(E_pot[:,1])**2
                tau[2] = np.where(C<=0.5*max(C))[0][0]
                err[2] = standdev(E_pot[:,1],tau[2],t_final)[0]
                mean[2] = standdev(E_pot[:,1],tau[2],t_final)[1]
                print("E, K, P are",mean,"error is",err)
      ########################
    
    elif state==0:
      ax.clear()
      ax.set_xlim3d([0.0, L])
      ax.set_xlabel('X')
      ax.set_ylim3d([0.0, L])
      ax.set_ylabel('Y')
      ax.set_zlim3d([0.0, L])
      ax.set_zlabel('Z')
      ax.set_title('Argon particles in a 3D box')
      
    #Fluctuations in the potential energy
    fluc_V=abs(initial_V-V)
    #print("V",V,"K",K,"E",total_E)

    # Checking conservation laws
    total_mom=check_mom(momenta)
    if sum(total_mom) > 1e-13:
        print("MOMENTA IS NOT CONSERVED")

    #if abs(initial_E-total_E) > fluc_V:
        #print("ENERGY IS NOT CONSERVED")

    return position, momenta, cnt, force, total_E, 

############################################################

 # Animation function
def animate(t_final, position, momenta, force, delta_t, initial_V, initial_E,cnt):

    position=iteration(position, momenta,
                       force, initial_E,
                       initial_V, delta_t, N, L,cnt)[0]
    ax.scatter(position[:,0], position[:,1], position[:,2], s=150)

    return position

############################################################

# Initial plot
def init_plot():

    position=init_position(A,L)
    ax.scatter(position[:,0], position[:,1], position[:,2], s=150, c='r')
    return position

############################################################
########################Main Run############################
############################################################

position = init_position(A, L)
momenta = init_momenta(T,N)
force, initial_E, initial_V, initial_K, initial_p = calc_force(position,momenta,N,Sig,Eps)
print("Initial energy", initial_E)
print("Initial potential", initial_V)

#############################################################################
#########################ANIMATION OR CORRELATION############################
#############################################################################

if state==0:
    cnt = 0.

    

    
    #Initializing the plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Setting the axes properties
    ax.set_xlim3d([0.0, L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, L])
    ax.set_zlabel('Z')

    ax.set_title('3D Test')

    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, animate, t_final,
                                       fargs=(position, momenta,
                                              force,delta_t, initial_V, initial_E,cnt),
                                       save_count=None,
                                       init_func=init_plot, blit=False)
    # Show the animation
    plt.show()

elif state==1:
    cnt = 0.
    #Make the iterative process run on its own
    for t in range(t_final):
      cnt += 1
      position, momenta, cnt = iteration(position,momenta,force,initial_E, initial_V,delta_t,N,L,cnt)[0 : 3]
      
        


#Print the global caltulation time
#print("Time(s):", time.clock() - start_time)
