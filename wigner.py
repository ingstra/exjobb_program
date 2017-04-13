from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from matplotlib import cm
from scipy.stats import norm
from scipy.misc import factorial
from scipy.special import eval_genlaguerre as L
from scipy.special import eval_laguerre as lag
import sympy.mpmath as mp
from mpl_toolkits.mplot3d import Axes3D

rho = np.loadtxt('rho.dat')

rho=np.array([[1, 0],[0,0]])
NFock = len(rho)

x=np.linspace(-3,3,500)
p=x

x, p = np.meshgrid(x,p)

laguerre = np.vectorize(L)

def W_mn(x,p,m,n):
    if m >= n:
        w = np.exp(-x*x-p*p)*np.power(-x+1j*p,m-n)*np.sqrt(np.power(2,m-n)*factorial(n)/factorial(m))* np.power(-1,m)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
    elif m < n: 
        w = np.exp(-x*x-p*p)*np.power(x+1j*p,n-m)*np.sqrt(np.power(2,n-m)*factorial(m)/factorial(n))* np.power(-1,n)*laguerre(m,n-m,2*x*x+2*p*p)/np.pi
    return w   


def wigner(rho,NFock,x,p):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return W

W1=wigner(rho,NFock,x,p)

def wigner2(rho,x,p):
    W=0
    for n in range(0,NFock):
        W=W+rho[n][n]*W_n(x,p,n)
        print rho[n][n],n
    return W

lag2=np.vectorize(lag)

def W_n(x,p,n):
    return 2*np.power(-1,n)*np.exp(-x*x-p*p)*lag(n,2*x*x+2*p*p)/np.pi

print rho



woho=wigner2(rho,x,p)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, p, np.real(W1),cmap=cm.coolwarm)

plt.show()


