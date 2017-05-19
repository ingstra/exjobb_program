from __future__ import division

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib as mlp
import numpy as np
from matplotlib import cm
from scipy.stats import norm
from scipy.misc import factorial
from scipy.special import eval_genlaguerre as L
from scipy import integrate

T=[0.2*i for i in range (1,26)]
Omega=[0.1*i for i in range (1,11)]


for n in T:
     for m in Omega:
          
          filename_im = "rho_im_T" + str(n) + "_Omega" + str(m)
          filename_re = "rho_re_T" + str(n) + "_Omega" + str(m)

          try:
               rho_real = np.loadtxt(filename_re)
               rho_im = np.loadtxt(filename_im)
          except IOError:
               print 'nofile: T:',n,'Omega:',m
               continue
          else:
               pass
