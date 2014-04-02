#! /usr/bin/env python

# This is some sandpit code for examining the results of Turing style reaction diffusion equations
# 

import warnings
import time
import sys
import numpy as np
import scipy.optimize
import scipy.integrate
import scipy.special
import random

import pdb

# This is just to remove some annoying convergence warnings that can really crowd the environment
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Here we define out reaction diffusion system
# t is time, y the state vector (corresponds to both Y_r and X_r in Turing's paper), and a the set of coefficients
def dfdt(xs, t, a, b, c, d, mu, nu):
    
    x = xs[1:len(xs)/2]
    y = xs[len(xs)/2:]

    # The vector xs consists of X_r in the first half, and Y_r in the second
    dxdt = np.zeros(len(x))
    dydt = np.zeros(len(y))
    
    for i in range(len(x)):
        # We take a modulo to take care of the circular arrangment of cells
        before = (i - 1) % len(x) / 2
        after = (i + 1) % len(x) / 2
        # These are equations (6.2) from the paper
        dxdt[i] += a * x[i] + b * y[i] + mu * (x[before] - 2 * x[i] + x[after])
        dydt[i] += c * x[i] + d * y[i] + nu * (x[before] - 2 * x[i] + x[after])

    result = np.append(dxdt, dydt) 
    return result

# Number of cells in the ring system
N = 20

# Now we run the code for a selected bunch of parameters.
# First try case a) from Section 8
a = 3 
b = 1
c = - 1
d = a
mu = 0.25
nu = mu 

# The initial conditions are perturbed just a little bit
xs_0 = np.ones(2*N) + 0.1 * np.random.random(2*N)

t = np.linspace(0., 10., 1000)

xs = scipy.integrate.odeint(dfdt, xs_0, t, args = (a,b,c,d,mu,nu))

print t
print xs


