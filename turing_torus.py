# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 15:40:57 2014

@author: hugobowne-anderson
"""

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
from pandas import *
from pylab import *
#%matplotlib inline
#%pylab inline

#set parameters for reaction-diffusion
I = 0.5
a = I-1 
b = 1
c = - 1
d = I
mu = 0.25
nu = mu
t = [0,0.01]
#t = range(0,100,1)

#define differential equations for reaction-diffusion
def dfdt(state, t, nbstate, a, b, c, d, mu, nu):
            dsdt = np.zeros(len(state))
            dsdt[0] = a*state[0] + b*state[1] + mu*(nbstate[0]-2*state[0])
            dsdt[1] = c*state[0] + d*state[1] + mu*(nbstate[1]-2*state[1])
            #for i in range(0,len(self.state)):
             #   dsdt[i] = 
            return dsdt
#integrate
            
#define class  Site: a "Site" is an individual cell in a Turing system, ok? i shuold probably have called it a cell. I 
# still need to include docstrings etc
class Site:
    def __init__(self,location,state):
        self.loc = location
        self.state = state
        
    def setNeighbours(self,neighbours):
        self.nb = neighbours
    
    def getNeighbourStates(self):
        self.nbstate = [0,0]
        for s in self.nb:
            for i in range(0,len(self.state)):
                self.nbstate[i] += s.state[i]
    
    def update(self):
        self.state = scipy.integrate.odeint(dfdt, self.state, t, args = (self.nbstate,a,b,c,d,mu,nu))[1]
    
    def __str__(self):
        return '<'  + str(self.loc) + '>'
        
#define class SiteCollection (yes, a colelction of Sites, good).
class SiteCollection:
    def __init__(self,sites):
        self.Sites = sites
        
    def update(self):
        for s in self.Sites:
            s.getNeighbourStates()
        for s in self.Sites:
            s.update()
            
    def __str__(self):
        return '<'  + str(self.Sites) + '>'
        
        
###############################################################################################################
#We the generate the data for a 2D torus of cells and save it to a CSV. I do not plot anything here
#because the best way to visualize it at the moment is as a dynamic heat map and I cannot seem to do that 
#in iPython. There is some R code in this folder to visualize it so make sure to check it out!
#I'll document this all further soon.
###############################################################################################################
"""
Created on Tue Jun 10 16:46:22 2014

@author: hugobowne-anderson
"""

#sys.setrecursionlimit(6000)

n = 100
L = [Site([i,j],[10,-10] + (n+8)*np.random.random(2)) for i in range(0,n) for j in range(0,n)]


for i in range(0,n**2):
    L[i].setNeighbours([ L[i/n*n + (i + 1)%n] , L[i/n*n + (i - 1)%n] , L[i-n] , L[(i+n)%(n**2)] ])
    
SC1 = SiteCollection(L)

s = 3000
stateX = np.zeros([n**2,s])
stateY = np.zeros([n**2,s])


import copy
for i in range(0,s):
    print(i)
    for j in range(0,n**2):
        stateX[j,i] = SC1.Sites[j].state[0]
        stateY[j,i] = SC1.Sites[j].state[1]
    SC1.update()


#  

#
#np.savetxt("stateY_n=150.csv", stateY, delimiter=",")
#
np.savetxt("stateX_n=150.csv", stateX, delimiter=",")
##        
