import integration as itg
import numpy as np

def Romberg(f,a,b,I):
    R= np.zeros([I,I])
    for i in range(I):
        m=0
        while m<=i:
            if m==0: #Ri1 using only trapezoid method
                N = 100 *2**i
                R[i,m] = itg.trapezoid(f,a,b,N) # doubling step size each row
            else:
                R_new = R[i,m-1]+1/(np.power(4,m-1))*(R[i,m-1]-R[i-1,m-1]) #correct error
                R[i,m]=R_new  
            m+=1
    return R


def f(x):
    return x**2

Romberg(f,0,2,3)