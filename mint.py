

# import needed libraries
from pylab import *
from scipy.integrate import dblquad

# define functions
# integrand
def sigma(y,x):            # NOTE ORDER OF ARGUMENTS! Required order for dblquad()
    return 1.3*exp(-x**2)  # gm cm^-3, doesn't actually depend on y

# lower limit of y integration
def ylower(x):             # just returns the constant value of 0
    return 0.

# upper limit of y integration
def yupper(x):
    return 2.*x/3.

# compute mass and display result
mass, error = dblquad(sigma,0.,3.,ylower,yupper)

print('The mass of the plate is %g grams' % mass)