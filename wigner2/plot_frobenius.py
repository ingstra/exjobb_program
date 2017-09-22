from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp


h =  np.loadtxt('frobenius1photon_dim2.dat')
#h =  np.loadtxt('frobenius.dat')


plt.semilogy(h[:,0],h[:,1],'k',linewidth=2)

plt.xlim([0,2000])

plt.xlabel(r'Iteration step')
plt.ylabel(r'$\Delta\rho$')
plt.tight_layout()
#plt.savefig('../figures/diff22.pdf',figsize=(10,10))
plt.show()
