#! /usr/bin/env python

#importing necessary libraries and a couple of things to get matplotlib working in this ipython notebook
import warnings
import time
import sys
import numpy as np
import scipy.optimize
import scipy.integrate
import scipy.special
import random

import matplotlib.pyplot as plt
from matplotlib import animation

from reaction_diffusion import *

#%matplotlib inline
#%pylab inline

a=-2
b=2.5
c=-1.25
d=1.5
mu=4
nu=2

#I=0.3
#a=I-1.
#b=1.
#c=-b
#d=I
#mu=1.
#nu=0

# Number of points
N = 100
# Length of ring
R = 1.0
# Final simulation time
T = 10.
# Number of time points
N_T = 101

# Define the Sites
L = [Site(i,[random.random(),random.random()], a, b, c, d, mu, nu) for i in range(0,N)]

for i in range(0,N):
    L[i].setNeighbours([L[i-1],L[(i+1)%N]])

collection = SiteCollection(L)
soln_t = np.linspace(0., T, N_T)
soln = collection.solve(soln_t)

#
# Here on is just animation code. Could be put in to a useful routine of some form...
#

xs = np.linspace(0., R, N)

fig = plt.figure(figsize=(8,8))
plt.xlim(0,R)
plt.ylim(soln.min(), soln.max())
plt.xlabel('x')
line1, = plt.plot([],[],'g-')
line2, = plt.plot([],[],'b-')
plt.legend([line1, line2], ["Morphogen 1", "Morphogen 2"])

def update(i, line1, line2):
    line1.set_data(xs, soln[0,i,:])
    line2.set_data(xs, soln[1,i,:])
    return line1, line2

# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, update, 
        frames=N, fargs=(line1, line2), interval=10)

# Save file as ScriptName_GitTag.mp4
import inspect, os, subprocess
exec_name =  os.path.splitext(os.path.basename(inspect.getfile(inspect.currentframe())))[0]
git_tag = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).replace('\n', '')

file_name = '{0}_{1}.mp4'.format(exec_name, git_tag)
print "Saving animation to", file_name

## UNCOMMENT THIS IF YOU WANT THE ANIMATION TO BE SAVED!
#anim.save(file_name, fps=24)
plt.show()
print 'what!'
