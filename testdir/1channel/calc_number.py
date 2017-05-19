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
import sympy.mpmath as mp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

rho = np.loadtxt('rho_re_.dat')

dim=4
x= np.arange(dim**2)
x = x.reshape((dim, dim))
N=np.zeros_like(x)
for i in range (0,dim):
    N[i,i]=i

mat = np.dot(N,rho)
trace = np.matrix.trace(mat)
print trace
    

