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

