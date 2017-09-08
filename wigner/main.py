import numpy as np
from numpy import kron
from scipy.linalg import expm
import time

start_time = time.time()
sigma_minus = np.matrix('0,0;1,0')
sigma_plus = np.matrix([[0,1],[0,0]])
rho0 = np.matrix([[1,0],[0,0]])
rho_vec0=np.transpose(rho0).A1

a=sigma_plus
adagger=sigma_minus

gamma=1
Omega=0.5
dt=1e-3
I = np.eye(2)
  
H =   -1j*gamma*np.matmul(sigma_plus,sigma_minus)/2 -1j*np.sqrt(gamma)*Omega*(sigma_plus-sigma_minus)
I = np.eye(2)

def L(c,H):
    return -1j*kron(I,H) + 1j*kron(np.conj(H),I) + kron(np.conj(c),c)

def M(c,lmbda):
    return np.sqrt(gamma)*( lmbda*kron(np.conj(c),I)-np.conj(lmbda)*kron(I,c) )

myfile = open('exact.dat', 'w')
rho=np.eye(2)

for i in range(1,10000):
    rho_vec = np.matmul(expm(i*dt*L(sigma_minus,H)),rho_vec0)
    rho = np.transpose(rho_vec.reshape(2,2))   
    result = np.trace(np.matmul(np.matmul(sigma_plus,sigma_minus),rho))
    myfile.write(str(i*dt) + '\t' + str(np.real(result)) + '\n')


print time.time() - start_time,' seconds to run.'
