from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm

#import sys
#trace = np.loadtxt('trace.dat')

traj = np.loadtxt('traj.dat')
exact = np.loadtxt('exact.dat')
current = np.loadtxt('current.dat')
#data = np.loadtxt('tst.dat')

dt = 1e-3
nruns = 1000

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)

#plt.hist(current,50)

#mufit, stdfit = norm.fit(current)
#print 'calculated mean: ',mu,'fitted mean: ',mufit
#print 'calculated std: ',sigma,'fitted std: ',stdfit


plt.plot(traj[:,0],traj[:,1])
plt.plot(exact[:,0],exact[:,1],'r--')



#x = mu + sigma*np.random.randn(1000)
#n, bins, patches = plt.hist(x, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--')


#plt.plot(trace[:,0],trace[:,1])


#plt.plot(traj[:,0],np.exp(-traj[:,0]),'r--')



def avg_Z(t):
    Omega = 1
    gamma = 1
    
    if 4*Omega < gamma:
        kappa = np.sqrt(np.power(gamma/4,2)-np.power(Omega,2))
    else:
        kappa = 1j*np.sqrt(-np.power(gamma/4,2)+np.power(Omega,2))
    
    K1 = np.power(Omega,2) / np.add(np.power(gamma,2) , 2*np.power(Omega,2) )
    hyp = np.cosh(kappa*t) + 3*gamma*np.sinh(kappa*t)/(4*kappa)

    return  K1*(1 - np.exp(-3*gamma*t/4)*hyp ) #- 1/2

def avg_Z_big_Omega(t):
    Omega = 5
    gamma = 1
    return -np.exp(-3*gamma*t/4)*np.cos(Omega*t)/2

def p_plus(t):
    Omega = 5
    gamma = 1
    return (1 - np.exp(-3*gamma*t/4)*np.cos(Omega*t))/2

#plt.plot(traj[:,0],avg_Z(traj[:,0]),'r--')

#plt.plot(data[:,0],p_plus(data[:,0]),'g')



plt.show()
