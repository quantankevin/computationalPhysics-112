import numpy as np
import matplotlib.pyplot as plt
import pickle

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
        self._tags = np.arange(N)
        self._time = 0.0
    pass
    def add_particles( self,masses, positions, velocities, accelerations):
        self._masses=np.concatenate((self._masses, masses), axis=0)
        self._positions=np.concatenate((self._positions, positions), axis=0)
        self._velocities=np.concatenate((self._velocities, velocities), axis=0)
        self._accelerations=np.concatenate((self._accelerations, accelerations), axis=0)
        self.nparticles+=np.shape(masses)[0]
    
    def plot_output(self, dim):
        fig = plt.figure()

        if dim == 2:
            ax = fig.add_subplot(111)
            ax.scatter(self.positions[:,0], self.positions[:,1])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            
        elif dim == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.positions[:,0], self.positions[:,1], self.positions[:,2])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            ax.set_zlabel('Z [code unit]')
        else:
            print("Invalid dimension!")
            return
    def save_data(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump({
                'masses': self.masses,
                'positions': self.positions,
                'velocities': self.velocities,
                'accelerations': self.accelerations
            }, f)      
    def load_data(self, filename):
        with open(filename, 'rb') as f:
            data = pickle.load(f)
            self.masses = data['masses']
            self.positions = data['positions']
            self.velocities = data['velocities']
            self.accelerations = data['accelerations']

    def save_data_txt(self, filename):
     with open(filename, 'w') as f:
        for i in range(self.nparticles):
            line = f"{self.masses[i][0]} {' '.join(map(str, self.positions[i]))} {' '.join(map(str, self.velocities[i]))} {' '.join(map(str, self.accelerations[i]))}\n"
            f.write(line)
    def load_data_txt(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
            self.nparticles = len(lines)
            self.masses = np.zeros((self.nparticles, 1))
            self.positions = np.zeros((self.nparticles, 3))
            self.velocities = np.zeros((self.nparticles, 3))
            self.accelerations = np.zeros((self.nparticles, 3))
        for i, line in enumerate(lines):
            data = line.split()
            self.masses[i] = float(data[0])
            self.positions[i] = list(map(float, data[1:4]))
            self.velocities[i] = list(map(float, data[4:7]))
            self.accelerations[i] = list(map(float, data[7:]))
    @property
    def time(self):
        return self._time
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
    def set_particles(self, pos, vel, acc):
        """
        Set particle properties for the N-body simulation

        :param pos: positions of particles
        :param vel: velocities of particles
        :param acc: accelerations of particles
        """
        self.positions = pos
        self.velocities = vel
        self.accelerations = acc
        return
    def output(self, filename):
        """
        Output particle properties to a file

        :param filename: output file name
        """
        masses = self.masses
        pos = self.positions
        vel = self.velocities
        acc = self.accelerations
        tags = self.tags
        time = self.time

        header = """
                ----------------------------------------------------
                Data from a 3D direct N-body simulation. 

                rows are i-particle; 
                coumns are :mass, tag, x ,y, z, vx, vy, vz, ax, ay, az

                NTHU, Computational Physics 

                ----------------------------------------------------
                """
        header += "Time = {}".format(time)
        np.savetxt(filename,(tags[:],masses[:,0],pos[:,0],pos[:,1],pos[:,2],
                            vel[:,0],vel[:,1],vel[:,2],
                            acc[:,0],acc[:,1],acc[:,2]),header=header)
        return
    


