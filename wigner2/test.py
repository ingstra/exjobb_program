import numpy as np
from qutip import *
import matplotlib as mlp
from matplotlib import cm
import matplotlib.pylab as plt

g = np.loadtxt('purity.dat')
plt.plot(g[:,0],g[:,1])

plt.show()
