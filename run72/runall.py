import subprocess
import numpy as np

T=[0.05*i for i in range (1,201)]
Omega=[0.05*i for i in range (1,21)]
print T

for n in T:
    for m in Omega:
        filename_im = "rho_im_T" + str(n) + "_Omega" + str(m)
        filename_re = "rho_re_T" + str(n) + "_Omega" + str(m)
        print str(n), str(m)
        output = subprocess.call(["./runprogram", str(n), str(m),filename_im,filename_re])
