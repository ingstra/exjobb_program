import numpy as np
from numpy import kron
from scipy.linalg import expm
import time

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

def L(c,H):
    I = np.eye(2)
    return -1j*kron(I,H) + 1j*kron(np.conj(H),I) + kron(np.conj(c),c)


myfile = open('g2.dat', 'w')

rho_vec = np.matmul(expm(1*dt*L(sigma_minus,H)),rho_vec0)

def evolve():
    for i in range(1,1000):
        rho_vec = np.matmul(expm(i*dt*L(sigma_minus,H)),rho_vec)
        rho = np.transpose(rho_vec.reshape(2,2))
        result = np.trace(np.matmul(np.matmul(sigma_plus,sigma_minus),rho))
        # myfile.write(str(i*dt) + '\t' + str(np.real(result)) + '\n')


P0 =kron(I,I)


def corr():
    rho_vec = np.matmul( expm(1*dt*L(a,H)), rho_vec0)
    rho = np.transpose(rho_vec.reshape(2,2))
    for i in range(2,20000):
        if (i*dt > 10):
            P = np.matmul( expm((i*dt-10.0)*L(a,H)), P0)
            g1 = np.real(np.trace(np.matmul(np.matmul(adagger,a),rho)))
            tmpvec = np.matmul(kron(np.conj(a),a),rho_vec)
            tmpvec = np.matmul(tmpvec,P)
            tmp= np.transpose(tmpvec.reshape(2,2))
            G2 = np.real(np.trace(np.matmul(adagger,a)*tmp))
            g2_norm = G2/g1**2
            myfile.write(str(i*dt-10) + '\t' + str(np.real(g2_norm)) + '\n')
    rho_vec = np.matmul( expm(i*dt*L(a,H)), rho_vec0)
    rho = np.transpose(rho_vec.reshape(2,2))

corr()


print time.time() - start_time,' seconds to run.'
