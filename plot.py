from __future__ import division

import matplotlib.pylab as plt
import numpy as np
import sys

# input file
filename = 'testfil'

data = np.loadtxt(filename)

#plt.hist(data)

plt.plot(data[:,0],data[:,1])

#plt.plot(data[:,0],np.exp(-2*data[:,0]),'y--')



def avg_Z(t):
    Omega = 5
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


plt.plot(data[:,0],avg_Z(data[:,0]),'r--')

#plt.plot(data[:,0],p_plus(data[:,0]),'g')



plt.show()
