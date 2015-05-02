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
import scipy.misc # This contains the image libraries
from PIL import Image # Or use this.... probably better than scipy.misc

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation

# This is a step-debugging library. Way useful.
import pdb

# This is a back-end that makes the code work for Mac OSX
#matplotlib.use('TkAgg')

# This is just to remove some annoying convergence warnings that can really crowd the environment
warnings.filterwarnings('ignore', category=RuntimeWarning)

def array_to_img(img_array, w, h):
    # reconstruct image here. Can't quite use reshape by itself...
    r = img_array[:len(img_array)/3].reshape(h,w)
    g = img_array[len(img_array)/3:2*len(img_array)/3].reshape(h,w)
    b = img_array[2*len(img_array)/3:].reshape(h,w)

    return np.dstack([r,g,b])

def laplacian(y, i, j, w, h, k):
    # Take Laplacian of vectorised 2-d field, with three species (i.e. k=0,1,2)
    return 0.25 * (y[(i-1)*w + j + k*w*h] + y[(i+1)*w + j + k*w*h] \
            + y[i*w + (j-1) + k*w*h] + y[i*w + (j+1) + k*w*h] - 4.*y[i*w + j + k*w*h]) 
    

# Here we define out reaction diffusion system
# t is time, y the state vector (corresponds to both Y_r and X_r in Turing's paper), and a the set of coefficients
def colour_diffuse(img, t, D, w, h):
    
    didt = np.zeros(img.shape)
    print "integrating at", t 
    for i in range(1,h-1):
        for j in range(1,w-1):
        
            ddx2 = np.array([laplacian(img, i, j, w, h, k) for k in range(3)])
            # This multiply gives us the diffusion forces as well as chemotactic forces
            forces = np.dot(D, ddx2)
            
            for k in range(3):
                # First we just add the diffusion part of the function, i.e. laplacian of each colour
                didt[i*w + j + k*w*h] = forces[k]
            
    return didt 

#img = misc.imread("james_hugo.jpg")
img = Image.open("james_hugo.jpg")
img_array = np.array(img)
# Flatten 3 2D channels in to one vector...
img_flat = np.append(np.append(img_array[:,:,0].flatten(), img_array[:,:,1].flatten()), img_array[:,:,2].flatten())

t = np.linspace(0., 1., 5)
dX = 500.0

# These are the diffusion terms
D_r = dX*0.01
D_g = dX*0.1
D_b = dX*0.02

# These are the chemotaxis terms: -ve means attraction, +ve means repulsion
n_rg = 0.      # Attraction of red towards green
n_rb = dX*0.05 # Attraction of red towards blue...
n_gr = -dX*0.03
n_gb = 0.
n_br = 0.
n_bg = 0.

D = dX * np.array([[D_r, n_rg, n_rb], [n_gr, D_g, n_gb], [n_br, n_bg, D_b]] )

# Need the width and the height of the image to process the vector
(h, w) = img_array[:,:,0].shape

img_states = scipy.integrate.odeint(colour_diffuse, img_flat, t, args = (D, w, h))
final_img_array = array_to_img(img_states[-1], w, h).astype(np.uint8)

final_img = Image.fromarray(final_img_array)
final_img.save('processed.jpg')
diff_img = Image.fromarray(final_img_array - img_array)
diff_img.save('difference.jpg')
