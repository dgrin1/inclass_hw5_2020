#!/usr/bin/env python

import numpy as np
from scipy.integrate import dblquad

#ran out of time for comments but you get the idea

def simpson_integration(f, a, b, N):
    h = (b-a)/N
    sum_odd = sum([f(a+k*h) for k in range(1, N, 2)])
    sum_even = sum([f(a+k*h) for k in range(2, N, 2)])
    return (1/3)*h*(f(a)+f(b)+4*sum_odd+2*sum_even)

def romberg(f, a, b, start_n, levels):
    Ilist = [ simpson_integration(f, a, b, n) 
              for n in [i*start_n for i in range(1, levels+1) ] ]
    arr = np.zeros((levels, levels))
    for i in range(levels):
        arr[i,0] = Ilist[i]
        for m in range(1, i+1):
            arr[i, m] = arr[i, m-1] + (
                1/(4**(m+1)-1))*(arr[i, m-1]-arr[i-1, m-1])
            
    return arr[-1, -1]


def sigma(y, x):
    return 1.3*np.exp(-x**2)


def center_of_mass():
    total_mass, error = dblquad(sigma, 0., 3., lambda x: 0, 
                                lambda x: (2.*x)/3.)


    xm, _ = dblquad(lambda y, x: x*sigma(y,x), 0, 3, 
                    lambda x: 0, lambda x: (2.*x)/3.)
    ym, _ = dblquad(lambda y, x: y*sigma(y,x), 0, 3, 
                    lambda x: 0, lambda x: (2.*x)/3.)

    xcm = (1/total_mass)*xm
    ycm = (1/total_mass)*ym

    print(xcm, ycm)

if __name__ == '__main__':
    f = lambda x: x**4-2*x+1 
    print(romberg(f, 0, 2, 20, 8))

    center_of_mass()
