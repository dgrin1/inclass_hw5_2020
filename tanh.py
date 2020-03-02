#Dan loves hg and detests github
from __future__ import division,print_function
from scipy import constants as sc
from numpy import loadtxt,array,dot,sqrt,linspace,log
from numpy import cos,sin,tanh
import numpy as np
from math import pi
import matplotlib.pyplot as plt

# Create area of x and tanh values (toy model of magnetization)
x=np.linspace(0,5,100,endpoint=True)
y=np.tanh(x)
h=0.01
#Simple forward difference derivative
def fdif_der(f,x,h):
#Compute forward difference
	s=f(x+h)-f(x-h)
#Divide by stepsize
	s/=(2*h)
	return s

# Enter interactive mode
#plt.ion()
#

plt.figure(1)
#Set up latex fonts
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','serif':['Times New Roman']})

#Multiple panels to visualize function, derivative, and error
plt.subplot(311)
plt.plot(x,y,'b',linewidth=3)
plt.ylim(0,1.3)
plt.ylabel('Magnetization'r'$M(T)$',fontsize=18)

#Plot Derivative
plt.subplot(312)
plt.plot(x,fdif_der(np.tanh,x,h),linewidth=1)
plt.plot(x,np.power(np.cosh(x),-2),'r*',markersize=10,linewidth=1)
plt.yscale('log')
plt.ylim([0,10])
plt.xlabel(r'$1/T~({\rm eV}^{-1})$',fontsize=18)
plt.ylabel(r'$\frac{dM}{dT}')

#Plot Fractional error
plt.subplot(313)
plt.plot(x,(np.abs((fdif_der(np.tanh,x,h)-np.power(np.cosh(x),-2))/np.power(np.cosh(x),-2))),'k',label='h')
plt.plot(x,(np.abs((fdif_der(np.tanh,x,h/10)-np.power(np.cosh(x),-2))/np.power(np.cosh(x),-2))),'r--',label='h/10')
plt.plot(x,(np.abs((fdif_der(np.tanh,x,h/100)-np.power(np.cosh(x),-2))/np.power(np.cosh(x),-2))),'b8',label='h/100')
plt.legend(loc=1)
plt.yscale('log')
plt.ylim([1.e-10,4.e-2])
plt.xlabel(r'$1/T~({\rm eV}^{-1})$',fontsize=18)

#Plot fractional error at a fixed x value, varying h
plt.figure(2)
Narr=np.logspace(1,3,100,endpoint='True')
plt.plot(Narr,np.abs((fdif_der(np.tanh,3,h/Narr)-np.power(np.cosh(3),-2))/np.power(np.cosh(3),-2)))
plt.xscale('log')
plt.yscale('log')
plt.ylim([1.e-20,1.e-2])
plt.show()








plt.show()


