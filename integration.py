import numpy as np
from gaussxw import gaussxw

def f(x):
	val=np.sin(x)
	return val

def riemann_dan(f,a,b,N):
#setting bin width
	delta_x=(b-a)/N
	s=0
#Loop over all y values
	for i in range(0,N,1):
#set x value at left edge of bin
		x=a+i*delta_x
#set function value of bin
		y=f(x)
#increment sum	
		s+=y
	s*=delta_x
	return s

def trapezoid(f,a,b,N):
	h = (b-a)/N
	s = 0

	for i in range (1,N):
		x = a + i*h
		s +=2*f(x)
	s += f(a)
	s += f(b)
	s *= h/2
	return s

def simpson(f,a,b,N):
	h = (b-a)/N
	s = f(a)
	loop = int(N/2)
	for i in range (0,loop):
		# start from x+1*h
		odd_i = 2*i+1
		x_odd = a + odd_i*h
		s += 4*f(x_odd)
		even_i = odd_i+1 #start from x+2*h
		if (even_i!=N): #not the last onw
			x_even = a + even_i*h
			s += 2*f(x_even)
	s += f(b)
	s *= h/3
	return s

def romberg(f,a,b,N):
    h = (b-a)/N
    s = 0
    return s


def gaussion(f,a,b,N):
    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    # rescale weight and x values
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    s = 0.0
    
    for k in range(N):
        s += wp[k]*f(xp[k])      
    return s