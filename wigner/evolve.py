import numpy as np
from numpy import kron
from numpy import sqrt
from scipy.linalg import expm
from scipy import integrate
import scipy as sp
import time
import functools
import itertools
import matplotlib.pylab as plt


start_time = time.time()
sigma_minus = np.matrix('0,0;1,0')
sigma_plus = np.matrix([[0,1],[0,0]])
rho0 = np.matrix([[0,0],[0,1]])
rho_vec0=np.transpose(rho0).A1
print rho0
gamma=1
Omega=0#.5
dt=1e-3

a=sigma_minus
adagger=sigma_plus
  
H =  -1j*np.sqrt(gamma)*Omega*(sigma_plus-sigma_minus) - 1j*gamma*np.matmul(sigma_plus,sigma_minus)/2

I = np.eye(2)

def L(c,H):
    return -1j*kron(I,H) + 1j*kron(np.conj(H),I) + kron(np.conj(c),c)

def M(c,l1,l2):
    return sqrt(gamma)*(l1+1j*l2)*kron(np.conj(c),I) - sqrt(gamma)*(l1 - 1j*l2)*kron(I,c)



myfile = open('exact.dat', 'w')

def evolve():
    for i in range(2,10000):
        rho_vec = np.matmul(expm(i*dt*L(sigma_minus,H)),rho_vec0)
        rho = np.transpose(rho_vec.reshape(2,2))
        result = np.trace(np.matmul(np.matmul(sigma_plus,sigma_minus),rho))
        myfile.write(str(i*dt) + '\t' + str(np.real(result)) + '\n')


P0 =kron(I,I)

corrfile = open('g2.dat', 'w')
t_start = 10

def corr():
    for i in range(1,20000):
        rho_vec = np.matmul( expm(i*dt*L(a,H)), rho_vec0)
        rho = np.transpose(rho_vec.reshape(2,2))
        if (i*dt > t_start):
            P = np.matmul( expm((i*dt-t_start)*L(a,H)), P0)
          
            g1 = np.real(np.trace(np.matmul(np.matmul(adagger,a),rho)))
            tmpvec = np.matmul(kron(np.conj(a),a),rho_vec)
           
            tmpvec = P*np.transpose(tmpvec)
           
            
            tmp= tmpvec.reshape(2,2)
            
            G2 = np.real(np.trace(np.matmul(np.matmul(adagger,a),tmp)))
            g2_norm = G2/g1**2
            
        
            corrfile.write(str(i*dt-10) + '\t' + str(np.real(g2_norm)) + '\n')
    

T=10           
I4 = np.eye(4) 
tau = 3 

def integrand(l1,l2,x,p):
    rhovec_t = np.matmul( expm(T*L(a,H)), rho_vec0)
    tmp = expm(tau*( L(a,H) - I4*(l1**2 + l2**2)/2 + M(a,l1,l2) ))
    tst = np.matmul(tmp,rhovec_t)
    tr_matrix = np.transpose(tst.reshape(2,2))
    tr = np.trace(tr_matrix)
    e=np.exp(2j*(p*l1 - x*l2))
    return tr*e/np.pi**2

def wigner(x,p):
    return integrate.dblquad(integrand, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf, args=(x,p))

x=np.linspace(-1,1,2)
p=x

X, P = np.meshgrid(x,p)

awesome_wigner=np.vectorize(wigner)
print wigner(0,0)
#W=awesome_wigner(X,P)
#print W
fig = plt.figure()
ax1 = fig.add_subplot(111)
#plt1 = ax1.contourf(x, p, W, 100, cmap=wcmap)


print time.time() - start_time,' seconds to run.'

