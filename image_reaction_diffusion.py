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

# Here we define out reaction diffusion system
# t is time, y the state vector (corresponds to both Y_r and X_r in Turing's paper), and a the set of coefficients
def colour_diffuse(img, t, D, w, h):
    
    # dfdt template is stupid - reconstruct img shape here
    #r = img[:len(img)/3].reshape(w,h)
    #g = img[len(img)/3:2*len(img)/3].reshape(w,h)
    #b = img[2*len(img)/3:].reshape(w,h)

    # The vector img consists of X_r in the first half, and Y_r in the second
    #drdt = np.zeros(r.shape)
    #dgdt = np.zeros(g.shape)
    #dbdt = np.zeros(b.shape)
    didt = np.zeros(img.shape)
    print "integrating at", t 
    for k in range(3):
        for i in range(1,h-1):
            for j in range(1,w-1):
                didt[i*w + j + k*w*h] = D[k] * 0.25 * (img[(i-1)*w + j + k*w*h] + img[(i+1)*w + j + k*w*h] \
                                        + img[i*w + (j-1) + k*w*h] + img[i*w + (j+1) + k*w*h] - 4.*img[i*w + j + k*w*h]) 

                #drdt[i,j] = D_r * 0.25 * (r[i-1,j] + r[i+1,j] + r[i,j-1] + r[i,j+1] - 4.*r[i,j]) 
                #dgdt[i,j] = D_g * 0.25 * (g[i-1,j] + g[i+1,j] + g[i,j-1] + g[i,j+1] - 4.*g[i,j]) 
                #dbdt[i,j] = D_b * 0.25 * (b[i-1,j] + b[i+1,j] + b[i,j-1] + b[i,j+1] - 4.*b[i,j]) 

    #result = np.append(np.append(drdt.flatten(), dgdt.flatten()), dbdt.flatten())
    return didt 

#img = misc.imread("james_hugo.jpg")
img = Image.open("james_hugo.jpg")
img_array = np.array(img)
# Flatten 3 2D channels in to one vector...
img_flat = np.append(np.append(img_array[:,:,0].flatten(), img_array[:,:,1].flatten()), img_array[:,:,2].flatten())

t = np.linspace(0., 1., 5)
dX = 0.05
D_r = dX*0.01
D_g = dX*0.1
D_b = dX*0.02

(h, w) = img_array[:,:,0].shape

img_states = scipy.integrate.odeint(colour_diffuse, img_flat, t, args = ([D_r, D_g, D_b], w, h))
final_img_array = array_to_img(img_states[-1], w, h).astype(np.uint8)
pdb.set_trace()
final_img = Image.fromarray(final_img_array)
final_img.save('processed.jpg')
