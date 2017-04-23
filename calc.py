import numpy as np
from scipy.integrate import quad
from scipy.misc import factorial
from scipy.special import eval_hermite as hermite


x=np.linspace(-7,11,3)

def psi1(n,x):
	return 1/(np.power(2,n)*factorial(n))*1/np.sqrt(2*np.pi)*np.exp(-np.power(x,2)/2)*np.power(hermite(n,x/np.sqrt(2)),2);


def psi2(n,x):
	return 1/(np.power(2,n)*factorial(n))*1/np.sqrt(np.pi)*np.exp(-np.power(x,2))*np.power(hermite(n,x),2);


n=2


#print np.trapz(psi1(5,x))
print 'psi1', quad(psi1,-np.inf,np.inf,args=n)

print 'psi2', quad(psi2,-np.inf,np.inf,args=n)

z=np.linspace(0,100,10)

print x
test = [ 1. ]
print psi1(n,test)

