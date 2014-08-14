import scipy.integrate
import numpy as np

import pdb

class Site(object):
    """ define class  Site: a "Site" is an individual cell in a Turing system, 
    ok? i shuold probably have called it a cell. I still need to include docstrings etc"""
    def __init__(self,location, state, a=0.5, b=1., c=-1., d=0.5, mu=0.25, nu=0.25):
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
        # Currently deprecated: this method isn't currently used
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
        print self.adjacency
        self.t = soln_times

    def build_adjacency_matrix(self):
        """ Takes the sites collection and builds an adjacency matrix """
        self.adjacency = np.zeros((len(self.Sites), len(self.Sites)), dtype=bool)
        
        # By no means is this the most efficient way of building the adjacency matrix
        # (e.g. we have symmetry...)
        for i in range(len(self.Sites)):

            for j in range(len(self.Sites)):
                if self.Sites[j] in self.Sites[i].nb:
                    self.adjacency[i,j] = True
    
    def build_param_vectors(self):
        self.a_vec = [site.a for site in self.Sites]
        self.b_vec = [site.b for site in self.Sites]
        self.c_vec = [site.c for site in self.Sites]
        self.d_vec = [site.d for site in self.Sites]
        self.mu_vec = [site.mu for site in self.Sites]
        self.nu_vec = [site.nu for site in self.Sites]

    def update(self):
        # Currently deprecated: this method isn't currently used
        for s in self.Sites:
            s.getNeighbourStates()
        for s in self.Sites:
            s.update()
    
    def solve(self):
        
        # Stupidly odeint can only operate on functions that work with 1-D arrays, so we have to use the flattened vector of site states.
        # The .transpose() function is used as the array comes out as a 2xN array, where as we want Nx2 for the flatten routine to have
        # the first state variable for the first half of the vector, and second for the second half.
        states = np.array([site.state for site in self.Sites]).transpose().flatten()
        self.build_param_vectors()
        
        self.solution = scipy.integrate.odeint(self.ensemble_dfdt, states, self.t)        
       
        for i in range(len(self.Sites)):
            # update the state information for each site, that is return it from the flattened solution vector
            self.Sites[i].state = np.array([self.solution[:,i], self.solution[:,i + len(self.Sites)]])

        return self.solution
         
    def ensemble_dfdt(self, states, t):
        dfdt = np.zeros(len(states))

        # Stupidly odeint can only operate on functions that work with 1-D arrays, so we have to use the flattened vector of site states 
        state_1 = states[:len(self.Sites)]
        state_2 = states[len(self.Sites):]
        
        dfdt[:len(self.Sites)] = self.a_vec * state_1 + self.b_vec * state_2 + self.mu_vec * (np.dot(self.adjacency, state_1) - self.adjacency.sum(1) * state_1)
        dfdt[len(self.Sites):] = self.c_vec * state_1 + self.d_vec * state_2 + self.nu_vec * (np.dot(self.adjacency, state_2) - self.adjacency.sum(1) * state_2)
        
        return dfdt
            
    def __str__(self):
        return '<'  + str(self.Sites) + '>'

