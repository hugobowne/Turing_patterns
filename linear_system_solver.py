# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 13:21:03 2014

@author: hugobowne-anderson
"""
#simple linear system modeling to form a basis for Turing patterning simulations: this code runs but is not finished

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.ion()


#linear system coefficients
a = 3
b = 1
c = -1
d = 3

# the system of ODES
def f(z,t):
    xi = z[0]
    yi = z[1]
    f0 = a*xi + b*yi
    f1 = c*xi + d*yi
    return [f0, f1]
    
#initial conditions
   
x0 = 5
y0 = 5
z0 = [x0,y0]
t  = np.linspace(0, 1, 1000)   # time grid

#solve the ODEs
soln = odeint(f,z0,t)
x = soln[:,0]
y = soln[:,1]

#plot that shit
plt.figure()
plt.plot(t,x, label = "morphogen 1")
plt.plot(t,y, label = "morphogen 2")
plt.xlabel('time (s)')
plt.ylabel('amount of morphogen')
