import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .particles import Particles
#  from numba import jit, njit, prange, set_num_threads

"""
The N-Body Simulator class is responsible for simulating the motion of N bodies



"""

class NBodySimulator:

    def __init__(self, particles: Particles):
        
        # TODO
        self.particles = particles
        self.time      = particles.time
        # set up the default simulation enviroments
        self.setup()
        return

    def setup(self, G=1,
                    rsoft=0.01,
                    method="RK4",
                    io_freq=10,
                    io_header="nbody",
                    io_screen=True,
                    visualization=False):
        """
        Customize the simulation enviroments.

        :param G: the graivtational constant
        :param rsoft: float, a soften length
        :param meothd: string, the numerical scheme
                       support "Euler", "RK2", and "RK4"

        :param io_freq: int, the frequency to outupt data.
                        io_freq <=0 for no output. 
        :param io_header: the output header
        :param io_screen: print message on screen or not.
        :param visualization: on the fly visualization or not. 
        """
        
        # TODO
        self.G = G
        self.rsoft = rsoft
        self.method = method
        if io_freq <= 0: io_freq = np.inf # no output
        self.io_freq = io_freq
        self.io_header = io_header
        self.io_screen = io_screen
        self.visualization = visualization
        return

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """

        # TODO
        self.dt = dt
        self.tmax = tmax

        nsteps = int(np.ceil(tmax/dt))
        time = self.time
        particles = self.particles


        # setup numerical meothd
        method = self.method
        if method.lower() == "euler":
            _advance_particles = self._advance_particles_Euler
        elif method.lower() == "rk2":
            _advance_particles = self._advance_particles_RK2
        elif method.lower() == "rk4":
            _advance_particles = self._advance_particles_RK4
        else:
            raise ValueError("Unknown method")

        # prepare for the output
        io_folder = "data_"+self.io_header
        Path(io_folder).mkdir(parents=True, exist_ok=True)

        # ===============================
        # Start the simulation
        # The main loop
        # ===============================
        for n in range(nsteps):

            # check if the time step exceeds the total time
            if (time+dt > tmax): dt = tmax - time

            # advance the particles
            particles = _advance_particles(dt, particles)

            # output data if needed
            if (n % self.io_freq == 0):
                if self.io_screen:
                    print("n=",n , "Time: ", time, " dt: ", dt)
                fn = self.io_header+"_"+str(n).zfill(6)+".dat"
                fn = io_folder+"/"+fn
                particles.output(fn)

            # update the time
            time += dt

        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        accelerations = np.zeros_like(positions)
        

        # TODO
        # Calculate the acceleration of the particles
        accelerations = np.zeros_like(positions)
        for i in range(nparticles):
            for j in range(nparticles):
                if i != j:
                    r = positions[j] - positions[i]
                    r_norm = np.linalg.norm(r)
                    a = self.G * masses[j] / (r_norm + self.rsoft)**3 * r
                    accelerations[i] += a

        # # return accelerations
        # accelerations = np.zeros_like(positions)
        # # invoke the kernel for acceleration calculation
        # accelerations = _calculate_acceleration_kernel(nparticles, masses, positions, accelerations, G, rsoft)

        return accelerations

    # @njit(parallel=True)
    # def _calculate_acceleration_kernel(nparticles, masses, positions, accelerations, G, rsoft):
    #     """
    #     Calculate the acceleration of the particles

    #     :param nparticles: int, the number of particles
    #     :param masses: ndarray, the masses of the particles
    #     :param positions: ndarray, the positions of the particles
    #     :param accelerations: ndarray, the accelerations of the particles
    #     :param G: float, the gravitational constant
    #     :param rsoft: float, the soften length
    #     """
    #     for i in prange(nparticles):
    #         for j in prange(nparticles):
    #             if i > j:
    #                 r = positions[j] - positions[i]
    #                 r_norm = np.linalg.norm(r)
    #                 a = G * masses[j] / (r_norm + rsoft)**3 * r
    #                 accelerations[i] += a
    #     return accelerations


    def _advance_particles_Euler(self, dt, particles):

        #TODO
        nparticles = particles.nparticles
        masses = particles.masses
        positions = particles.positions
        velocities = particles.velocities
        accelerations = self._calculate_acceleration(nparticles,masses, positions)

        # do the Euler update
        positions += velocities*dt
        velocities += accelerations*dt

        # update the particles
        particles.set_particles(positions, velocities, accelerations)





        return particles

    def _advance_particles_RK2(self, dt, particles):

        # TODO
        nparticles = particles.nparticles
        mass = particles.masses

        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # do the RK2 update
        pos2 = pos + vel*dt
        vel2 = vel + acc*dt
        acc2 = self._calculate_acceleration(nparticles, mass, pos2)

        pos2 = pos2 + vel2*dt
        vel2 = vel2 + acc2*dt

        # average
        pos = 0.5*(pos + pos2)
        vel = 0.5*(vel + vel2)
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # update the particles
        particles.set_particles(pos, vel, acc)
        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        #TODO
        nparticles = particles.nparticles
        mass = particles.masses

        # y0
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)
        
        dt2 = dt/2
        # y1
        pos1 = pos + vel*dt2
        vel1 = vel + acc*dt2
        acc1 = self._calculate_acceleration(nparticles, mass, pos1)

        # y2
        pos2 = pos + vel1*dt2
        vel2 = vel + acc1*dt2
        acc2 = self._calculate_acceleration(nparticles, mass, pos2)

        # y3
        pos3 = pos + vel2*dt
        vel3 = vel + acc2*dt
        acc3 = self._calculate_acceleration(nparticles, mass, pos3)

        # do the RK4 update
        pos = pos + (vel + 2*vel1 + 2*vel2 + vel3)*dt/6
        vel = vel + (acc + 2*acc1 + 2*acc2 + acc3)*dt/6
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # update the particles
        particles.set_particles(pos, vel, acc)
        return particles



if __name__ == "__main__":
    
    pass