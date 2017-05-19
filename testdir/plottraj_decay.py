from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

#import sys
#trace = np.loadtxt('trace.dat')

traj = np.loadtxt('traj1.dat')
exact = np.loadtxt('exact.dat')


dt = 1e-3
nruns = 10000

gamma = 1

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)


plt.plot(traj[:,0],traj[:,1],'b',linewidth=2)
plt.plot(exact[:,0],exact[:,1],'r--',linewidth=2)

#plt.plot(traj[:,0],np.exp(-gamma*traj[:,0]),'r:',linewidth=2)

#plt.plot(trace[:,0],trace[:,1])


plt.xlabel(r'$t$')
plt.ylabel(r'$\braket{\sigma_+\sigma_-}$')
plt.tight_layout()
plt.savefig('../figures/ntrajs1000_decay.pdf',figsize=(10,10))

plt.show()
