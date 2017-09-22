from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
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
from qutip import *

from wigner_cmap import w_cmap

rho_real = np.loadtxt('rho_re_.dat')
rho_im = np.loadtxt('rho_im.dat')

rho2 = rho_real + 1j*rho_im

rho=np.array([[0, 0],[0,1]])
#rho2 = np.array([[0.0004702, 0.008834],[0.008834,0.9995]])
#rho2 = np.array([[0.9888,-0.3629E-09],[-0.3629E-09,0.1122E-01]])

NFock = len(rho)

x=np.linspace(-3,3,500)
p=x

X, P = np.meshgrid(x,p)
dx=x[1]-x[0]
dp=dx
laguerre = np.vectorize(L)

def W_mn(x,p,m,n):
    if n>=m:
          w = np.exp(-x*x-p*p)* np.power(-1,m)*np.sqrt(np.power(2,n-m)*factorial(m)/factorial(n))*np.power(x-1j*p,n-m)*laguerre(m,n-m,2*x*x+2*p*p)/np.pi
    elif n < m: 
        w = np.exp(-x*x-p*p)* np.power(-1,n)*np.sqrt(np.power(2,m-n)*factorial(n)/factorial(m))*np.power(x+1j*p,m-n)*laguerre(n,m-n,2*x*x+2*p*p)/np.pi
    return w   

def wigner_calc(x,p,rho):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.real(W)


def wigner_calc_abs(x,p,rho):
    W=0
    for m in xrange(0,NFock):
        for n in xrange(0,NFock):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.abs(np.real(W))

def wignerint(x,p,rho):
    return wigner_calc_abs(x,p,rho) - wigner_calc(x,p,rho)

print rho

W=wigner_calc(X,P,rho)
W2 = wigner_calc(X,P,rho2)

fig = plt.figure()

plt.xlabel(r'$x$')
plt.ylabel(r'$W(x)$')
#plt.title(r'$N=1$ Fock state')

range= [-6,6]
#print integrate.nquad(wignerint,[range,range])
#print integrate.nquad(wigner_calc,[range,range],)
#print integrate.dblquad(wignerint,-np.inf, np.inf,lambda x: -np.inf, lambda x: np.inf)
Wabs=np.abs(W)
int1 = np.trapz(Wabs-W)*dp
int2=np.trapz(int1)*dx
print int2

plt.tight_layout()
plt.plot(x,W2[:,250],'b',label='Reconstruction',linewidth=2)
plt.plot(x,W[:,250],'r--',label='Exact',linewidth=2)

plt.legend()
print p[249:251]

#xvec = np.linspace(-5,5,200)
#W = wigner(rhoq, x, p)

#ax1 = fig.add_subplot(111)


wcmap = w_cmap(W,shift=0)

#plt1 = ax1.contourf(x, p, W, 100, cmap=wcmap)
#cb1 = fig.colorbar(plt1, ax=ax1)

#plt.savefig('../figures/wigner_1photon_comp.pdf',figsize=(10,10))


plt.show()
