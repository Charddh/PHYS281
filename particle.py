"""
Created on Fri Nov 18 11:08:23 2022

@author: c.d.hamilton@lancaster.ac.uk

contains class Particles and methods to update properties of the particles.
"""

import numpy as np
from astropy.constants import G

class Particles:
    """
    Class to set up default properties of particles to be included in simulation
    
    Parameters:
    -------------
    position, velocity, acceleration: numpy arrays
        Arrays containing the 3 dimensional positions, velocities and accelerations of the particles.
    name: string
        Name of the particle
    mass: float
        Mass of the particle
    KE, PE, LM, AM: floats
        Current kinetic energies, gravitational potential energies, linear momenta and angular momenta for the particles
    """
    
    def __init__(
        self,
        position=np.array([0, 0, 0], dtype=float),
        velocity=np.array([0, 0, 0], dtype=float),
        acceleration=np.array([0, 0, 0], dtype=float),
        name='Ball',
        mass=1.0,
        KE = 0,
        PE = 0,
        LM = 0,
        AM = 0,
        LA = np.array([0, 0, 0], dtype=float)):
        self.name = name
        self.position = np.array(position,dtype=float)
        self.velocity = np.array(velocity,dtype=float)
        self.acceleration = np.array(acceleration,dtype=float)
        self.mass = mass
        self.KE = 1
        self.PE = 1
        self.LM = 1
        self.AM = 1
        self.LA = np.array([0, 0, 0], dtype=float)
        
    def __str__(self):
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
            self.name, self.mass,self.position, self.velocity, self.acceleration
        )
           
    def update(self, deltaT, method):   
        """
        Method to update the positions and velocities of the particles
        
        Parameters:
        -------------
        deltaT: float
            The time interval for the simulation
        method: integer
            Used to choose between the three methods:
                1: Euler's method
                2: Euler-Cromer method
                3: Verlet algorithm
        """        
        if method == 1:
            self.position = np.array(self.position + (deltaT * self.velocity))
            self.velocity = np.array(self.velocity + (deltaT * np.array(self.acceleration)))        
        if method == 2:       
            self.velocity = np.array(self.velocity + (deltaT * np.array(self.acceleration)))
            self.position = np.array(self.position + (deltaT * self.velocity))
        if method == 3:
            self.position = np.array(self.position + (np.array(self.velocity * deltaT)) + np.array((0.5 * self.LA * (deltaT ** 2))))
            self.velocity = np.array(self.velocity + np.array((0.5 * (self.acceleration + self.LA)) * deltaT))
    
    def updateGravitationalAcceleration(self, bodies):
        """
        Method to update the accelerations of the particles
        
        Parameters:
        -------------
        bodies: list
            List containing the instances of the Particles class
        """   
        GravitationalAcceleration = 0
        for i in bodies:
            if i != self:
                accmag = -(G * i.mass) / np.square((np.linalg.norm(self.position - i.position)))
                accdirection = (self.position - i.position) / (np.linalg.norm(self.position - i.position))
                acc = accmag * accdirection
                GravitationalAcceleration += acc
        self.LA = self.acceleration
        self.acceleration = GravitationalAcceleration + 0.
        
    def updateKE(self):
        """
        Method to update the kinetic energies of the particles
        """
        self.KE = 0.5 * self.mass * np.linalg.norm(self.velocity)**2 + 0.
    
    def updatePE(self, bodies):
        """
        Method to update the gravitational potential energies of the particles
        
        Parameters:
        -------------
        bodies: list
            List containing the instances of the Particles class
        """
        PotentialE = 0
        for i in bodies:
            if i != self: 
                PotentialE += -(G.value * i.mass * self.mass) / (np.linalg.norm(self.position - i.position))
        self.PE = PotentialE + 0.
    
    def updateLM(self):
        """
        Method to update the linear momenta of the particles
        """
        self.LM = self.mass * self.velocity + 0.
        
    def updateAM(self):
        """
        Method to update the angular momenta of the particles
        """
        self.AM = np.cross(self.position, self.LM)