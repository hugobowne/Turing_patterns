import scipy.integrate
import numpy as np

import pdb

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

class Site(object):
    """ define class  Site: a "Site" is an individual cell in a Turing system, 
    ok? i shuold probably have called it a cell. I still need to include docstrings etc"""
    def __init__(self,location,state, a=0.5, b=1., c=-1., d=0.5, mu=0.25, nu=0.25):
        self.loc = location
        self.state = state
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.mu = mu
        self.nu = nu
        
    def setNeighbours(self, neighbours):
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

class SubSite(Site):
    """ Here I define a class SubSite in order to recursively divide Sites into smaller Sites and to respect
    neighbour-neighbour interactions """
    def __init__(self,location,state,parent):
        self.loc = location
        self.parent = parent

class SiteCollection(object):
    """ define class SiteCollection (yes, a colelction of Sites, good). """
    def __init__(self,sites, soln_times=np.linspace(0.0,1.0,101)):
        self.Sites = sites
        self.build_adjacency_matrix()

    def build_adjacency_matrix(self):
        """ Takes the sites collection and builds an adjacency matrix """
        self.adjacency = np.zeros((len(self.Sites), len(self.Sites)), dtype=bool)
        
        # By no means is this the most efficient way of building the adjacency matrix
        # (e.g. we have symmetry...)
        for i in range(len(self.Sites)):

            for j in range(len(self.Sites)):
                if self.Sites[j] in self.Sites[i].nb:
                    self.adjacency[i,j] = True
    
    #def get_site_vector(self):
    #    return [site.State for site in self.Sites]

    def update(self):
        for s in self.Sites:
            s.getNeighbourStates()
        for s in self.Sites:
            s.update()
    
    def better_update(self):
        self.state = scipy.integrate.odeint(dfdt, self.state, t, args = (self.nbstate,a,b,c,d,mu,nu))[1]        
         
    def ensemble_dfdt(self):
        dsdt = np.zeros(len(state))
        
        a_vec = [site.a for site in self.Sites]
        b_vec = [site.b for site in self.Sites]
        c_vec = [site.c for site in self.Sites]
        d_vec = [site.d for site in self.Sites]
        mu_vec = [site.mu for site in self.Sites]
        nu_vec = [site.nu for site in self.Sites]
        
        states = [site.State for site in self.Sites]
        
        pdb.set_trace()
        dfdt[0] = a_vec * states[0,:] + b * states[1,:] + mu * (dot(self.adjacency * states[0,:]) - self.adjacency.sum(1) * states[0,:])
        dfdt[1] = c_vec * states[0,:] + b * states[1,:] + mu * (dot(self.adjacency * states[1,:]) - self.adjacency.sum(1) * states[1,:])

        
        #for i in range(0,len(self.state)):
         #   dsdt[i] = 
        return dfdt
        
            
    def __str__(self):
        return '<'  + str(self.Sites) + '>'

