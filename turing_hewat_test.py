# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 16:08:54 2014

@author: hugobowne-anderson
"""

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
#define class  Site: a "Site" is an individual cell in a Turing system, ok?
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