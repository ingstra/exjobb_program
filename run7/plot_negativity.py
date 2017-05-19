from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib as mlp
import numpy as np
from matplotlib import cm
from scipy.stats import norm
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from wigner_cmap import w_cmap


T=[0.2*i for i in range (1,101)]
Omega=[0.1*i for i in range (1,26)]


N=np.loadtxt('negmatrix.dat')

print len(T),len(Omega),N.shape

fig = plt.figure()
ax1 = fig.add_subplot(111)
wcmap = w_cmap(N,shift=0)
plt1 = ax1.contourf(Omega,T, N/2, 100, cmap=cm.PuRd)
plt.xlim([0.1,1])
plt.ylim([0.1,8.5])

cb1 = fig.colorbar(plt1, ax=ax1)


plt.show()
