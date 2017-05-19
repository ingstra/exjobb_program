import numpy as np
#from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt
from scipy.special import eval_genlaguerre as L
from scipy.misc import factorial

from visualize import matrix_histogram




rho_real = np.loadtxt('rho_re_.dat')
rho_im = np.loadtxt('rho_im.dat')

rho = rho_real + 1j*rho_im

print(rho)

N = len(rho)

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

def wigner_calc(x,p):
    W=0
    for m in xrange(0,N):
        for n in xrange(0,N):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.real(W)


def wigner_calc_abs(x,p):
    W=0
    for m in xrange(0,N):
        for n in xrange(0,N):
            W = W + rho[m][n]*W_mn(x,p,m,n)     
            
    return np.abs(np.real(W))

def wignerint(x,p):
    return wigner_calc_abs(x,p) - wigner_calc(x,p)

W=wigner_calc(X,P)




                                   
xlabels = [r'$\bra{0}$',r'$\bra{1}$',r'$\bra{2}$',r'$\bra{3}$',r'$\bra{4}$',r'$\bra{5}$'] 

ylabels = [r'$\ket{0}$',r'$\ket{1}$',r'$\ket{2}$',r'$\ket{3}$',r'$\ket{4}$',r'$\ket{5}$']
fig, ax = matrix_histogram(rho[0:4,0:4], xlabels, ylabels,limits=[-1,1])
ax.view_init(azim=-55, elev=25)

#rho_ss = steadystate(rho, [np.sqrt(0.1) * a, np.sqrt(0.4) * b.dag()])



#fig, ax = hinton(rho, xlabels=ylabels, ylabels=xlabels)



#plt.savefig('../figures/rho_1channels_T5_Omega2_2.pdf',figsize=(10,10))

plt.show()
