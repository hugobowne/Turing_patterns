#! /usr/bin/env python

#importing necessary libraries and a couple of things to get matplotlib working in this ipython notebook
import copy
import warnings
import time
import sys
import numpy as np
import scipy.optimize
import scipy.integrate
import scipy.special
import random

#from pandas import *
from pylab import *
import matplotlib.pyplot as plt
from reaction_diffusion import *

#%matplotlib inline
#%pylab inline

a=0.5
b=1.
c=-1.
d=0.5
mu=0.25
nu=0.25

I = 0.3
a=I-1.
b=1.
c=-b
d=I
mu=1.
nu=0

## Define a collection of 20 Sites such that they form an icosahedron (brute-force method!)
#define the Sites
L = [Site(i,[random(),random()], a, b, c, d, mu, nu) for i in range(0,20)]

T = 1.

#defien the neighbours according to an icosahedral net i will include in an accompanying .tex
L[0].setNeighbours([L[1],L[16],L[4]])
L[1].setNeighbours([L[0],L[2],L[18]])
L[2].setNeighbours([L[3],L[1],L[5]])
L[3].setNeighbours([L[2],L[7],L[19]])
L[4].setNeighbours([L[5],L[0],L[8]])
L[5].setNeighbours([L[4],L[6],L[2]])
L[6].setNeighbours([L[7],L[5],L[9]])
L[7].setNeighbours([L[6],L[1],L[3]])
L[8].setNeighbours([L[9],L[4],L[12]])
L[9].setNeighbours([L[8],L[10],L[6]])
L[10].setNeighbours([L[11],L[9],L[13]])
L[11].setNeighbours([L[10],L[15],L[7]])
L[12].setNeighbours([L[13],L[6],L[16]])
L[13].setNeighbours([L[12],L[14],L[10]])
L[14].setNeighbours([L[15],L[13],L[17]])
L[15].setNeighbours([L[14],L[19],L[11]])
L[16].setNeighbours([L[17],L[12],L[0]])
L[17].setNeighbours([L[16],L[18],L[14]])
L[18].setNeighbours([L[19],L[17],L[1]])
L[19].setNeighbours([L[18],L[13],L[15]])

# For now I'm ignoring the subsite business, just going to test this one on 
collection = SiteCollection(L, np.linspace(0., T, 101))
soln = collection.solve()

plt.subplot(1,2,1)
for site in collection.Sites:
    plt.plot(collection.t, site.state[0], collection.t, site.state[1])

plt.subplot(1,2,2)
plt.plot(range(len(collection.Sites)), soln[-1,:len(collection.Sites)])

plt.show()
