import numpy as np
from numpy import kron
from numpy import sqrt
from scipy.linalg import expm
from scipy import integrate
import scipy as sp
import time
import functools
import itertools
from matplotlib import cm
import matplotlib.pylab as plt
from wigner_cmap import w_cmap

start_time = time.time()
sigma_minus = np.matrix('0,0;1,0')
sigma_plus = np.matrix([[0,1],[0,0]])
rho0 = np.matrix([[0,0],[0,1]])
rho_vec0=np.transpose(rho0).A1
gamma=1
Omega=0#.5
dt=1e-3

a=sigma_minus
adagger=sigma_plus
  
H =   - 1j*gamma*np.matmul(sigma_plus,sigma_minus)/2 -1j*np.sqrt(gamma)*Omega*(sigma_plus-sigma_minus)

I = np.eye(2)

def L(c,H):
    return -1j*kron(I,H) + 1j*kron(np.conj(H),I) + kron(np.conj(c),c)

def M(c,l1,l2):
    return sqrt(gamma)*(l1+1j*l2)*kron(np.conj(c),I) - sqrt(gamma)*(l1 - 1j*l2)*kron(I,c)

T=  0.00001     
I4 = np.eye(4) 
tau = 3 

def integrand(l1,l2,x,p):
    rhovec_t = np.matmul( expm(T*L(a,H)), rho_vec0)
    tmp = expm(tau*( L(a,H) - I4*(l1**2 + l2**2)/2 + M(a,l1,l2) ))
    tst = np.matmul(tmp,rhovec_t)
    tr_matrix = np.transpose(tst.reshape(2,2))
    tr = np.trace(tr_matrix)
    e=np.exp(2j*(p*l1 - x*l2))
    return np.real(tr*e/np.pi**2)

def wigner(x,p):
    return integrate.dblquad(integrand, -5, 5, lambda x: -5, lambda x:5, args=(x,p))

lim=6
def wigner2(x,p):
    return integrate.nquad(integrand,[[-lim,lim],[-lim,lim]],args=(x,p))[0]

n=6
x=np.linspace(-lim,lim,n)
p=x


X, P = np.meshgrid(x,p)

awesome_wigner=np.vectorize(wigner2)

#print wigner(0,0)

#test= wigner2(0,0)
#print test


#for i in range(0,n-1):
#    for j in range(0,n-1):
#        W[i,j]=wigner(x[i],p[j])


W=awesome_wigner(X,P)
print W
fig = plt.figure()
ax1 = fig.add_subplot(111)
#wcmap = w_cmap(W,shift=0)
plt1 = ax1.contourf(x, p, W, 100, cmap=cm.coolwarm)


print time.time() - start_time,' seconds to run.'

plt.show()
