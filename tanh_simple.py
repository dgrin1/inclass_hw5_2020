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
plt.ion()

#plt.figure(1)
#Set up latex fonts
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','serif':['Times New Roman']})

#Multiple panels to visualize function, derivative, and error
plt.subplot(211)
plt.plot(x,y,'b',linewidth=3)
plt.ylim(0,1.3)
plt.ylabel('Magnetization'r'$M(T)$',fontsize=18)

#Plot Derivative
plt.subplot(212)
plt.plot(x,fdif_der(np.tanh,x,h),linewidth=1)
plt.plot(x,np.power(np.cosh(x),-2),'r*',markersize=10,linewidth=1)
plt.yscale('log')
plt.ylim([0,10])
plt.xlabel(r'$1/T~({\rm eV}^{-1})$',fontsize=18)
plt.ylabel(r'$\frac{dM}{dT}')
plt.show()


