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

start_time = time.time()
sigma_minus = np.matrix('0,0;1,0')
sigma_plus = np.matrix([[0,1],[0,0]])
rho0 = np.matrix([[1,0],[0,0]])
rho_vec0=np.transpose(rho0).A1

gamma=1
Omega=0.5
dt=1e-3

a=sigma_minus
adagger=sigma_plus
  
H =  -1j*np.sqrt(gamma)*Omega*(sigma_plus-sigma_minus) - 1j*gamma*np.matmul(sigma_plus,sigma_minus)/2

I = np.eye(2)
rho=np.array([[0, 0],[0,1]])

def integrand(l1,l2,x,p):
    tr = np.trace(np.matmul(rho,expm((l1+1j*l2)*adagger-(l1-1j*l2)*a)))
    e=np.exp(2j*(p*l1 - x*l2))
    return tr*e/np.pi**2

lim=100

def wigner(x,p):
    return integrate.nquad(integrand,[[-lim,lim],[-lim,lim]],args=(x,p))
  #  return integrate.dblquad(integrand, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf, args=(x,p))



N=1000

def fftwigner(j,k,x,p):
    y=0
    for l in range(-N/2,N/2):
        for n in range(-N/2,N/2):
           y=y+ np.exp(4j*pi*k*l/N)*np.exp(-4j*pi*j*n/N)*np.trace(np.matmul(rho,expm((l+1j*n)*adagger - (l-1j*n)*a) ) )/pi**2
    return y

def fftwigner2(j,k,x,p):
    y=0
    for l in range(0,N-1):
        for n in range(0,N-1):
           y=y+ np.exp(4j*pi*k*l/N)*np.exp(-4j*pi*j*n/N)*np.trace(np.matmul(rho,expm((l+1j*n)*adagger - (l-1j*n)*a) ) )/(N**2 * pi**2)
    return y


M=61
dt=0.2
dx=0.1
def f(l1,l2):
    return np.cos(sqrt(l1**2+l2**2)/2)


def W_kl(k,l):
    y=0
    for m in range(0,M-1):
        for n in range(0,M-1):
            y = y + f(m*dt,n*dt)*np.exp(2*pi*1j*k*m/M)*np.exp(-2*pi*1j*n*l/M)
    return y/(4*pi**2)

def W_kl2(k,l):
    y=0
    for m in range(-M/2,M/2):
        for n in range(-M/2,M/2):
            y = y + np.exp(4*pi*1j*k*m/M)*np.exp(-4*pi*1j*n*l/M)*np.trace(np.matmul(rho,expm((m+1j*n)*adagger - (m-1j*n)*a)))
    return y/pi**2


print W_kl(0,0)
print W_kl2(0,0)
print fftwigner(0,0,0,0)
print fftwigner2(0,0,0,0)

x=np.linspace(-30,30,M)
print x
p=x
X, P = np.meshgrid(x,p)
W=W_kl(X*dx,P*dx)





#l1=1
#l2=2
#print np.cos(np.sqrt(l1**2 + l2**2))
#print np.trace(np.matmul(rho,expm((l1+1j*l2)*adagger - (l1-1j*l2)*a) ))

#awesome_wigner=np.vectorize(wigner)

#W=awesome_wigner(X,P)
#print W
#wcmap = w_cmap(W,shift=0)
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
#surf = ax1.plot_surface(X*dx, P*dx, np.real(W),cmap=cm.coolwarm)


print time.time() - start_time,' seconds to run.'

#plt.show()
