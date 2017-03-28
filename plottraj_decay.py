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
current = np.loadtxt('current.dat')

dt = 1e-3
nruns = 10000

gamma = 1

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)


plt.plot(traj[:,0],traj[:,1],'k',linewidth=2)
plt.plot(exact[:,0],exact[:,1],'g--',linewidth=2)

plt.plot(traj[:,0],np.exp(-gamma*traj[:,0]),'r:',linewidth=2)

#plt.plot(trace[:,0],trace[:,1])


plt.tight_layout()
plt.savefig('testplot',figsize=(20,10))

plt.show()
