import numpy as np

def relax(func, x_init):

    x = x_init
    while func(x)!=x:
        x=func(x)

    return x

def func1(x):

    return np.exp(1-x**2)

y = relax(func1, .9)
print(y)
