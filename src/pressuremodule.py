import numpy as np
import matplotlib.pyplot as plt


def pres(lpnum,cutoff,n,prestime,toten,kenarray,potarray,Ttarg,density,N):
  averagedp = np.zeros((lpnum-cutoff)/n, dtype=float)
  for p in range((lpnum-cutoff)/n):
    averagedp[p] = np.mean(prestime[n*p:n*(p+1)])
     
  np.savetxt("presdataT%srho%sN%sn%slpnum%sprestime.txt"%(Ttarg,density,N,n,lpnum),prestime)
  np.savetxt("presdataT%srho%sN%sn%slpnum%saveragedp.txt"%(Ttarg,density,N,n,lpnum),averagedp)
  np.savetxt("totenT%srho%sN%sn%slpnum%s.txt"%(Ttarg,density,N,n,lpnum),toten)
  np.savetxt("kenT%srho%sN%sn%slpnum%s.txt"%(Ttarg,density,N,n,lpnum),kenarray)
  np.savetxt("potT%srho%sN%sn%slpnum%s.txt"%(Ttarg,density,N,n,lpnum),potarray)



  #fig1 = plt.plot(kenarray)
  #fig2 = plt.plot(toten)
  #fig3 = plt.plot(toten-kenarray)
  fig4 = plt.plot(averagedp)
  plt.show()

  mnp = np.mean(prestime)
  sdomp = np.std(prestime)/np.sqrt(len(prestime))
  mnavp = sum(averagedp)/len(averagedp)
  sdomavp = np.std(averagedp)/np.sqrt(len(averagedp))

  print mnp, sdomp, mnavp, sdomavp
  return
