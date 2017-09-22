import numpy as np
from numpy import kron
from numpy import sqrt,pi
from scipy.linalg import expm
from scipy import integrate
import scipy as sp
import time
import functools
import itertools
import matplotlib.pylab as plt
from wigner_cmap import w_cmap
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from wigner_cmap import w_cmap

start_time = time.time()
sigma_minus = np.matrix('0,0;1,0')
sigma_plus = np.matrix([[0,1],[0,0]])
rho0 = np.matrix([[1,0],[0,0]])
rho_vec0=np.transpose(rho0).A1


a=sigma_minus
adagger=sigma_plus
  

I = np.eye(2)
rho=np.array([[0, 0],[0,1]])

n=100
x=np.linspace(-3,3,n)
p=x
dx=x[1]-x[0]
dp=dx

X, P = np.meshgrid(x,p)
tau=3

def wigner(x,p):
    return (2/(pi*tau))*np.trace(expm(-2/tau*np.absolute(tau*sigma_minus-x*I-1j*p*I)**2)*rho)

awesome_wigner=np.vectorize(wigner)
W=awesome_wigner(X,P)

lim=6
#print integrate.nquad(awesome_wigner,[lim,lim]) 

int1 = np.trapz(W)*dp
int2=np.trapz(int1)*dx
print int2
print wigner(0,0)

#wcmap = w_cmap(W,shift=0)
fig = plt.figure()
#ax1 = fig.add_subplot(111, projection='3d')
ax1 = fig.add_subplot(111)
#surf = ax1.plot_surface(X, P, np.real(W),cmap=cm.coolwarm)
wcmap = w_cmap(W,shift=0)
plt1 = ax1.contourf(x, p, W, 100, cmap=wcmap)
cb1 = fig.colorbar(plt1, ax=ax1)
print time.time() - start_time,' seconds to run.'

plt.show()
