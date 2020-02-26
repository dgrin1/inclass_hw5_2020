#from __future__ import print_function,division
import numpy as np



#Trapezoidal rule , defined as a sum, note that the function 
#f can be passed, making this a versatile routine that can be used with any function
def trape(a,b,f,N): #python lets you pass whole functions!!!!!
	h=(b-a)/N
	s=(f(a)+f(b))*h
	for i in range(1,N):
		s+=2.*f(a+i*h)*h
	s/=2.0
	return s

#Simpson's rule	
def simp(f,a,b,n): #python lets you pass whole functions
#initialize integral as zero
	int=0
#check that number of points is even
	if n%2 != 0:	
		print('Uneven number of points')
		h=0
		int=0
	else:
		h=(b-a)/float(n) # set step size
		int=(f(a)+f(b))*h/3.0e0 #outer and inner boundaries
		for i in range(1,n+1,2):
			int+=f(a+float(i)*h)*4.e0*h/3.e0 #add odd points
		for i in range(2,n,2):
			int+=f(a+float(i)*h)*2.e0*h/3.e0 #add even points
	return int,h

#Simpsons rule with an error estimate
def simpreal(f,a,b,n): #python lets you pass whole functions
	i1,h1=simp(f,a,b,n) # estimate integral with simpson's rule
	i2,h2=simp(f,a,b,2*n) #estimate again with half the step size
	err=(i2-i1)/15.e0 #use this to estimate the error
	i=i2
#	print(i1,i2)
	return i2,err,h1,h2	 #output integral better version, error, and two stepsizes
#t,err,stepa,stepb=simpreal(f,a,b,N)
	


#Compute integrand using trigonometric transform
def arg(kappa,theta):
	out=(kappa+np.tan(theta))**2.
	out*=np.exp(-np.tan(theta))
	out/=np.exp(kappa)-np.exp(-np.tan(theta))
	out*=np.power(np.cos(theta),-2.e0)
	return out

def frac_simpreal(kappa,N):
#define integrand for integral with binding energy
	def integrand(theta):
		return arg(kappa,theta)
#define integral for probability normalization		
	def integrand_botto(theta):
		return arg(0.,theta)
#use simpsons rule
	u,eu,a,b=simpreal(integrand,1.e-8,np.pi/2,N)
#Halve stepsize to get an error estimate
	v,ev,a,b=simpreal(integrand_botto,1.e-8,np.pi/2,N)

	return u/v,np.abs(u/v-(u+eu)/(v-ev))









#Make plot off error as a function of N
#Generate arrays
errarr=[]
errarr_g=[]
kappa=13.6/10.
nr=np.logspace(1,5,20)
nr=nr.astype(int)

for i in range(nr.size):
	if nr[i]%2==1:	
		nr[i]+=1
		
for i in nr:
	s,t=frac_simpreal(13.6/8.,i)
	errarr.append(t)


 #Make plot
import matplotlib.pyplot as plt
plt.ion()
plt.subplot(211)

plt.plot(nr,errarr)
plt.xscale('log')
plt.yscale('log')

#Arrays and graphics for second plot plots ionization fraction as a function of temperature
karr=np.arange(1,1000.6,0.05)
s,e=frac_simpreal(13.6/karr,70)
plt.subplot(212)
plt.ylim([1.e-4,2])
plt.xlim([1.e-0,1.e3])
plt.plot(karr,s)
plt.xscale('log')
plt.yscale('log')
plt.show()


















