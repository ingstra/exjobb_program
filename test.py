import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

rho= fock_dm(2)
print(rho)

xvec = np.linspace(-5,5,200)

W = wigner(rho, xvec, xvec)


def plot_wigner_2d_3d(psi):
    #fig, axes = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12, 6))
    fig = plt.figure(figsize=(17, 8))
    
    ax = fig.add_subplot(1, 2, 1)
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=6);

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    plot_wigner(psi, fig=fig, ax=ax, projection='3d', alpha_max=6);
    
    
    return fig

plot_wigner_2d_3d(rho)




plt.show()
