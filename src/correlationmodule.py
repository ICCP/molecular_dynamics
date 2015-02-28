import numpy as np
import matplotlib.pyplot as plt
import pylab

def cor(npdim,N,distances,nbins,finalbins,plotflag,Ttarg,density):
#find the correlation function of the system
  dmax = np.sqrt(5.)*npdim
  bin_vec = np.zeros(nbins, dtype = float)

  for i in range(N):
    for j in range(i+1,N):
      bin_num = int(distances[i][j]*nbins/dmax)
      if bin_num < nbins/2.:
        bin_vec[bin_num] = bin_vec[bin_num] + 1
      
      


#normalize based on the volume of the radial shell
#keeping order R^2*dR -> V=4*PI*R^2*dR
  dR = dmax/nbins
  Rin = np.zeros((nbins), dtype=float)
  for bin_num in range(nbins):
    Rout = (bin_num + 2)*dR
    Rin[bin_num] = (bin_num + 1)*dR
    bin_vec[bin_num] = bin_vec[bin_num]/(4*np.pi*Rin[bin_num]**2*dR)
    
  
  finalbins[0] = 0.0001 # to make sure that the plot starts at r=0  
  if plotflag == 1:
    plt.plot(Rin,finalbins)#,width=dR)
    plt.xlim([0,max(Rin)*(5/8.)])
    plt.ylabel('g(r)')
    plt.xlabel('r')
    plt.title('Correlationfunction')
    plt.text(max(Rin)*(1/4.),max(finalbins)-0.1,r'$T$=%s, $\rho$=%s'%(Ttarg,density)) # x and y values for position of text are choosen for 864 particles
    plt.show()
    print max(Rin)
    plt.savefig('correlationfunctionT%sRho%s.jpg'%(Ttarg,density))
    np.savetxt("cordataT%srho%sN%sRin.txt"%(Ttarg,density,N),Rin)
    np.savetxt("cordataT%srho%sN%sfinbin.txt"%(Ttarg,density,N),finalbins)
    
    
 
  return bin_vec
  

