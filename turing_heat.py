# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 15:46:41 2014

@author: hugobowne-anderson
"""

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
from pandas import *

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation

# This is a step-debugging library. Way useful.
import pdb

# This is a back-end that makes the code work for Mac OSX
#matplotlib.use('TkAgg')

# This is just to remove some annoying convergence warnings that can really crowd the environment
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Here we define out reaction diffusion system
# t is time, y the state vector (corresponds to both Y_r and X_r in Turing's paper), and a the set of coefficients
def dfdt(xs, t, a, b, c, d, mu, nu):
    
    x = xs[:len(xs)/2]
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

# The initial conditions are perturbed just a little bit
xs_0 = np.ones(2*N) + 0.1 * np.random.random(2*N)


n = 1200 # number of time steps
tt = 12 # total time of simulation
# This currently animates the solution to (c) in Section 8. Blue line is X_r, green line is Y_r
# Play with parameters to see other behaviour
t = np.linspace(0., tt, n)
I = 0.5
a = I-1 
b = 1
c = - 1
d = I
mu = 1
nu = 0 

#integrate
xs = scipy.integrate.odeint(dfdt, xs_0, t, args = (a,b,c,d,mu,nu))

#set up initial frame
y1 = xs[1, :2*N]
y2 = np.reshape(y1,(2,N))
#quad = plt.pcolormesh(y2,vmin = np.min(xs), vmax = np.max(xs),shading='gouraud')
quad = plt.pcolormesh(y2,vmin = np.min(xs), vmax = np.max(xs))
plt.colorbar()


#show initial frame: top row contains concentrations of morphogen X in the cells, bottom row of morphogen Y
plt.ion()
plt.show()

#now we iterate over time and show each frame <= this is a MOVING IMAGE, MAN
for j in t[1:]:
    i = (n-1)/tt*j-1 #indexing time with integers
    y1 = xs[i, :2*N] # vector of morphogen concentrations
    y2 = np.reshape(y1,(2,N))
    plt.title('time: %.2f'%j)
    quad.set_array(y2.ravel())
    plt.draw()
    
plt.ioff()
plt.show()
