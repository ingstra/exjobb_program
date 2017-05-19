import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

g = np.loadtxt('g2_Omega0_5.dat')
g2 = np.loadtxt('g2_Omega2.dat')

plt.plot(g[:,0],g[:,1],linewidth=2,label='$\Omega=0.5$')
plt.plot(g2[:,0],g2[:,1],'r--',linewidth=2,label='$\Omega=2$')
plt.ylim([0,2])
plt.xlim([0,10])


params = {'legend.fontsize': 20}
plt.rcParams.update(params)
plt.legend()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$g^{(2)}(\tau)$')
plt.tight_layout()
plt.savefig('../figures/g2.pdf',figsize=(10,10))

plt.show()
