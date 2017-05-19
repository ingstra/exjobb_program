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
import matplotlib.ticker as ticker
from wigner_cmap import w_cmap
from matplotlib import rc


T=[0.2*i for i in range (1,41)]
Omega=[0.1*i for i in range (1,11)]


N=np.loadtxt('negmatrix.dat')

N2=N[0:41,0:11]/2

print len(T),len(Omega),N2.shape

fig = plt.figure()
ax1 = fig.add_subplot(111)
wcmap = w_cmap(N2,shift=0)
#plt1 = ax1.contourf(Omega,T, N/2, 100, cmap=cm.PuRd)

imgplot=plt.imshow(N2,aspect='auto',interpolation='none', origin='lower')

#plt.xlim([0.1,1])
#plt.ylim([0.1,8.5])

#cb1 = fig.colorbar(plt1, ax=ax1)

imgplot.set_cmap('Purples')

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.1))
ax1.xaxis.set_major_formatter(ticks)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.2))
ax1.yaxis.set_major_formatter(ticks)
plt.colorbar()
plt.xlabel('$\Omega$')
plt.ylabel('$T$')
plt.tight_layout()
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.savefig('../figures/map.pdf',figsize=(10,10))
plt.show()
