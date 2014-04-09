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
N = 50

# The initial conditions are perturbed just a little bit
xs_0 = np.ones(2*N) + 0.1 * np.random.random(2*N)

# This currently animates the solution to (c) in Section 8. Blue line is X_r, green line is Y_r
# Play with parameters to see other behaviour
t = np.linspace(0., 10., 500)
I = 0.5
a = I-1 
b = 1
c = - 1
d = I
mu = 0.25
nu = mu 


xs = scipy.integrate.odeint(dfdt, xs_0, t, args = (a,b,c,d,mu,nu))

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-np.max(abs(xs)), np.max(abs(xs))))
line1, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 1, N)
    y1 = xs[i, :N]
    y2 = xs[i, N:]
    line1.set_data(x, y1)
    line2.set_data(x, y2)
    return line1, line2

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(t), interval=1, blit=False)
plt.show()
