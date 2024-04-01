import numpy as np
import matplotlib.pyplot as plt


class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N:int=100):
        self.nparticles = N
        self._masses = np.ones((N,1))
        self._positions = np.zeros((N,3))
        self._velocities = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))
    pass
    def add_particles( self,masses, positions, velocities, accelerations):
        self._masses=np.concatenate((self._masses, masses), axis=0)
        self._positions=np.concatenate((self._positions, positions), axis=0)
        self._velocities=np.concatenate((self._velocities, velocities), axis=0)
        self._accelerations=np.concatenate((self._accelerations, accelerations), axis=0)
        self.nparticles+=np.shape(masses)[0]
    
    def output(self, dim):
        if dim ==2:
            plt.scatter(self._positions[:,0],self._positions[:,1])
        elif dim ==3:
            plt.scatter(self._positions[:,0],self._positions[:,1],self.positions[:,2])
    pass

    @property
    def tags(self):
        return self._tags
    @tags.setter
    def tags(self, some_tags):
        if len(some_tags) != self.nparticles:
            print("Number of particles does not match!")
            raise ValueError    
        self._tags = some_tags
        return
   
    @property
    def masses(self):
        return self._masses
    @masses.setter
    def masses(self, some_masses):
        if np.shape(some_masses) !=  np.shape(self._masses):
 
            print("Number of particles does not match!")
            raise ValueError    
        self._masses = some_masses
        return
    
    @property
    def velocities(self):
        return self._velocities
    @velocities.setter
    def velocities(self, some_velocities):
        if np.shape(some_velocities) != np.shape(self._velocities):
            print("Number of particles does not match!")
            raise ValueError    
        self._velocities = some_velocities
        return
    
    @property
    def accelerations(self):
        return self._accelerations
    @accelerations.setter
    def accelerations(self, some_accelerations):
        if np.shape(some_accelerations) != np.shape(self._accelerations):
            print("Number of particles does not match!")
            raise ValueError    
        self._accelerations= some_accelerations
        return
    
    @property
    def positions(self):
        return self._positions
    @positions.setter
    def positions(self, some_positions):
        if np.shape(some_positions) != np.shape(self._positions):
            print("Number of particles does not match!")
            raise ValueError    
        self._positions= some_positions
        return


