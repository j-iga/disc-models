"""
A class for stellar objects, useful for keeping track of parameters during modelling.

Author: Jordan Iga, The University of Edinburgh
        Imperial College London UROP Student
"""

import numpy as np

# CONSTANTS https://physics.nist.gov/cuu/Constants/
C = 29979245800  # Speed of light in a vacuum (cm/s)
H_P = 6.62607015E-27  # Planck constant in cm^2 g/s
K_B = 1.380649E-16  # Boltzmann constant in cm^2 g/s^2/K
G = 6.67430E-8  # Newton's constant of Gravitation in cm^3/s^2
M_SOLAR = 1.989E33  # Solar Mass in grams
PARSEC = 3.0857E18  # Parsec in cm
R_SOLAR = 6.957E10  # Solar radius in cm


class Star(object):
    """
    A class to define a stellar object.

    Properties:
    label: A label for the star.
    mass: The Mass of the object in in Solar masses.
    radius: The radius of the object in Solar radii.
    t_eff: The effective temperature of the object in K.
    distance: The distance of the system from the sun in parsec.

    Methods:
    _init_
    set_radius: Determines and sets the radius of a object from its surface gravity.
    radius_cgs: Returns the radius of the star in cm.
    mass_cgs: Returns the mass of the star in grams.
    distance_cgs: Returns the distance of the star in cm.

    """

    def __init__(self, label, mass, radius, t_eff, distance):

        self.label = label
        self.mass = mass
        self.radius = radius
        self.t_eff = t_eff
        self.distance = distance

    def set_radius(self, log_g):

        self.radius = np.sqrt((G*self.mass)/(10**log_g))

    def radius_cgs(self):

        return self.radius*R_SOLAR

    def mass_cgs(self):

        return self.mass*M_SOLAR

    def distance_cgs(self):

        return self.distance*PARSEC



