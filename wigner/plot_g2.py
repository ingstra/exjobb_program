from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from scipy.misc import factorial

exact = np.loadtxt('g2.dat')



plt.plot(exact[:,0],exact[:,1],linewidth=2)


plt.tight_layout()

plt.show()
