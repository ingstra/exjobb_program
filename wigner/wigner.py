from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib as mlp
import numpy as np
from matplotlib import cm
from scipy.stats import norm
from scipy.misc import factorial
from scipy.special import eval_genlaguerre as L
from scipy.special import eval_laguerre as lag
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from qutip import *

from wigner_cmap import w_cmap

lim=3
n=15
x=np.linspace(-lim,lim,n)
y=x

X,Y=np.meshgrid(x,y)

file = np.loadtxt('wigner.dat')

W=file
#print W
wcmap = w_cmap(W,shift=0)

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt1 = ax1.contourf(x, y, W, 100, cmap=wcmap)
cb1 = fig.colorbar(plt1, ax=ax1)


plt.show()