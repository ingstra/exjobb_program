from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial
import sympy.mpmath as mp

current = np.loadtxt('current.dat')

dt = 1e-3
nruns = 1300

gamma = 2

mu = 0
variance = dt*nruns
sigma = np.sqrt(variance)

values, bins, _ = plt.hist(current,20,normed=1,histtype='bar',edgecolor = "black")
area = sum(np.diff(bins)*values)
#print bins

mufit, stdfit = norm.fit(current)
#print 'fitted mean: ',mufit #, 'calculated mean: '
print 'fitted var: ',stdfit**2 ,'calculated var: ',variance


range = 3
x= np.linspace(-range,range,100)

n=0
h_poly = np.frompyfunc(mp.hermite,2,1)

const = 1/(np.sqrt(2**n * factorial(n)) * (np.pi)**0.25 )
rest = const* np.exp(-x**2/2)
tot  = (rest*h_poly(n,x))**2
plt.plot(x,tot,'r--',linewidth=2)
plt.xlabel(r'Integrated current')

plt.tight_layout()
#print np.trapz(tot,x)
#sigma = 1
#s = np.random.normal(0, sigma, 1000)
#plt.hist(s,bins=50,normed=1)




plt.savefig('../figures/photocurrent_vacuum.pdf',figsize=(10,10))



plt.show()
