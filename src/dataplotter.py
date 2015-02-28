import numpy as np
import matplotlib.pyplot as plt

# Importing the data

Rin = np.loadtxt("cordataT0.5rho1.2N864Rin.txt")
finalbins = np.loadtxt("cordataT0.5rho1.2N864finbin.txt")
prestime = np.loadtxt("presdataT3rho0.3N864n100lpnum10000prestime.txt")
averagedp = np.loadtxt("presdataT3rho0.3N864n100lpnum10000averagedp.txt")
kenarray = np.loadtxt("kenT3rho0.5N864n100lpnum10000.txt")
potarray = np.loadtxt("potT1rho0.88N864n100lpnum10000.txt")
totenarray = np.loadtxt("totenT1rho0.88N864n100lpnum10000.txt")


# Values for the different parameters, need to be equal to those mentioned in the filenames.

Ttarg = 1
density = 0.88
n = 100 #width of pressure block
N = 864 #number of particles
lpnum = 10000
cutoff = int(lpnum/5.)


# Plotmodules

def corplotter(Rin,finalbins,Ttarg,density):
  plt.plot(Rin,finalbins)#,width=dR)
  plt.xlim([0,6])
  plt.ylabel('g(r)')
  plt.xlabel('r')
  plt.title('Correlationfunction')
  plt.text(4.2,max(finalbins)-0.1,r'$T$=%s, $\rho$=%s'%(Ttarg,density), fontsize=18) # x and y values for position of text are choosen for 864 particles
  plt.show()
  return
  
  
def presplotter(prestime,averagedp,n):
  fig = plt.figure()
  ax1 = fig.add_subplot(2,1,1)
  ax2 = fig.add_subplot(1,1,1)
  fig.subplots_adjust(hspace=.35)
  ax1.plot(range(len(prestime)),prestime)
  ax2.plot(range(len(averagedp)),averagedp)
  ax1.set_ylabel('pressure')
  ax1.set_xlabel('time')
  ax2.set_ylabel('pressure')
  ax2.set_xlabel('time')
  ax1.set_title('pressure evolution')
  ax2.set_title('pressure evolution with pressure blocks')
  ax2.text(1,2.715,'blocksize = %s'%(n),fontsize=12)
  plt.show()
  return
  
def energyplotter(kenarray,potarray,totenarray):
  plt.plot(range(len(kenarray)),kenarray)
  plt.plot(range(len(potarray)),potarray)
  plt.plot(range(len(totenarray)),totenarray)
  plt.title('Energy evaluation')  
  plt.xlabel('time')
  plt.ylabel('energy')
  plt.text(6800,-6500,r'N=%s, T=%s, $\rho$=%s'%(N,Ttarg,density), fontsize = 12)
  plt.show()
  return
  
def renormplotter(kenarray,totenarray):
  plt.plot(range(len(kenarray[0:2500])),kenarray[0:2500])
  plt.plot(range(len(totenarray[0:2500])),totenarray[0:2500])
  plt.title('Energy renormalization')  
  plt.xlabel('time')
  plt.ylabel('energy')
  plt.text(1600,-5000,r'N=%s, T=%s, $\rho$=%s'%(N,Ttarg,density), fontsize = 12)
  plt.show()
  return
 
 
# Calling the plot you want to make  

#renormplotter(kenarray,totenarray)
#energyplotter(kenarray,potarray,totenarray)
#presplotter(prestime,averagedp,n)



# Calculating the errors of the pressure and temperature

mnp = np.mean(prestime)
sdp = np.std(prestime)
sdpb = np.std(averagedp)
sdompb =  np.std(averagedp)/np.sqrt(len(averagedp))
mnT = np.mean(kenarray[cutoff:9999])*(2/(3.*N))
sdT = np.std(kenarray[cutoff:9999])*(2/(3.*N))
sdomT = sdT/np.sqrt(len(kenarray[cutoff:9999]))

print 'mnT=',mnT,'sdT=',sdT, 'sdomT=',sdomT, 'mnp=',mnp,'sdp=',sdp,'sdpb=',sdpb,'sdompb=',sdompb 
