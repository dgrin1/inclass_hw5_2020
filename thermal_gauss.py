from gaussxw import *
import numpy as np

#Compute integrand using trigonometric transform
def arg(kappa,theta):
	out=(kappa+np.tan(theta))**2.
	out*=np.exp(-np.tan(theta))
	out/=np.exp(kappa)-np.exp(-np.tan(theta))
	out*=np.power(np.cos(theta),-2.e0)
	return out

#This time with Gaussian integration
def frac_gaussreal(kappa,N):
#define integrand for integral with binding energy
	def integrand(theta):
		return arg(kappa,theta)
#define integrand for probability normalization		
	def integrand_botto(theta):
		return arg(0.,theta)
#import x coordinates and weights
	xgrid,wgrid=gaussxwab(N,1.e-8,np.pi/2)
	
	u=sum(wgrid*integrand(xgrid))
	v=sum(wgrid*integrand_botto(xgrid))
#Do twice as many evaluations	
	N*=2
	xgrid,wgrid=gaussxwab(N,1.e-8,np.pi/2)
	un=sum(wgrid*integrand(xgrid))
	vn=sum(wgrid*integrand_botto(xgrid))
#estimate errors
	eu=np.abs(un-u)
	ev=np.abs(vn-v)
	return u/v,np.abs(u/v-(u+eu)/(v-ev))
	
#test
s,t=frac_gaussreal(13.6/8.0,100)
print(s)
print(t)