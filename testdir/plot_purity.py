import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

g = np.loadtxt('purity.dat')
plt.ylim([0,1.1])
plt.plot(g[:,0],g[:,1],'-')


#plt.savefig('../figures/1channels_T5_Omega2_2.pdf',figsize=(10,10))
plt.show()
