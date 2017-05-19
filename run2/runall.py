import subprocess
import numpy as np

T=[0.2*i for i in range (1,26)]
Omega=[0.1*i for i in range (1,11)]

for n in T:
    for m in Omega:
        filename_im = "rho_im_T" + str(n) + "_Omega" + str(m)
        filename_re = "rho_re_T" + str(n) + "_Omega" + str(m)
        print str(n), str(m)
        output = subprocess.call(["./run2", str(n), str(m),filename_im,filename_re])
