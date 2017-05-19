import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

g = np.loadtxt('g2_Omega5.dat')
plt.plot(g[:,0],g[:,1])
plt.ylim([0,2])
plt.xlim([0,10])
plt.show()
