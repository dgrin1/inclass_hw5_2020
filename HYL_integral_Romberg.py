# 02/20/2019   
# inclass work for integration -- Romberg Integral
import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Trapezoid Sum (using equation in chapter 5 p142 of Newman)
def TrapezoidSum(f,a,b,N):
    delta=float(abs(b-a)/N)
    sum=0.0
    for i in range(1,N):
        # taking care of the case if a>b
        if b > a:
            sum+=f(a+i*delta)
        else:
            sum+=f(b+i*delta)
    sum=sum+(1/2)*(f(a)+f(b))
    sum*=delta
    # taking care of the case if a>b
    if b > a:
        return sum
    else:
        return sum*(-1.0)

# Romberg integation implemented according to chapter 5.4 p160 of Newman
def RombergIntegral(f,a,b,N):
    # Define a dictionary for storing the estimated integral R_(i,m)
    R={} 
    # Define the first two estimates of the integral using trapezoidal rule
    R[1,1]=TrapezoidSum(f,a,b,N)
    R[2,1]=TrapezoidSum(f,a,b,2*N)
    # Define a variable to update the iteration number for TrapezoidSum()
    Nmin=2*N
    # Define variables for keeping track of the indexes of the estimated Ri,m
    i=2
    m=1
    epsilon=1e-16 # target accuracy
    error=1e0 # variable for checking the computational error
    while abs(error) > epsilon: 
        # Conditional case that calculates R_i,m whenever i is larger than m
        if m<i:    
            # Defind the variable R_i,m+1 for iteration
            R[i,m+1]=R[i,m]+((1.0)/((4**m)-1))*(R[i,m]-R[i-1,m])
            # Compute the error using eqn-5.49 from chapter 5.4 p161 of Newman
            error = ((1.0)/((4**m)-1))*(R[i,m]-R[i-1,m])
            m+=1
        # Conditional case that calculates R_i,m whenever m is equal to i and i is upgraded 1 higher
        else:
            # Upgrade the value of the parameter i and reset the value of the parameter m
            i+=1
            m=1
            # Update the minimum iteration number for TrapezoidSum()
            Nmin*=2
            # Update the elements in the dictionary R_(i,m)
            R[i,m]=TrapezoidSum(f,a,b,Nmin)
    return R[i,m],error,i,Nmin/N

# test
print("Uisng Trapezoid Integral with 10000 iterations gives the result:", TrapezoidSum(lambda x: x**4-2*x+1,0,2,10000))
print("Using Romberg Integral gives the result:", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[0],"with a trapezoidal integration level of", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[2],"and an iteration number of",RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[3])
print("The error using Romberg Integral is:", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[1])

        
