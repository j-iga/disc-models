"""

Contains functions useful for the modelling of spectra from debris discs of White Dwarfs (WDs).
This program uses the flat disc model seen in M. Jura 2003.
The WhiteDwarf class is used.

Author: Jordan Iga, The University of Edinburgh
        Imperial College London UROP Student

"""


import numpy as np
from scipy import integrate


# CONSTANTS
C = 29979245800  # Speed of light in a vacuum (cm/s)
H_P = 6.62607015E-27  # Planck constant in cm^2 g/s
K_B = 1.380649E-16  # Boltzmann constant in cm^2 g/s^2/K
G = 6.67430E-8  # Newton's constant of Gravitation in cm^3/s^2
M_SOLAR = 1.989E33  # Solar Mass in grams
PARSEC = 3.0857E18  # Parsec in cm
R_SOLAR = 6.957E10  # Solar radius in cm
AU = 1.49598073E13  # Astronomical Unit in cm


def t_func(r, wd):
    """
    The radial temperature function for a flat, optically thick disc. (1) from M. Jura 2003.

    :param r: The distance of the annulus from the WD at the centre of the ring
    :param wd: A WhiteDwarf object
    :return t: The temperature of the annulus at r.
    """
    r = r*R_SOLAR
    t = (2/(3*np.pi))**(1/4) * (wd.radius/r)**(3/4) * wd.t_eff

    return t


def flux(freq, t_in, t_out, inc, wd):
    """
    Determines the flux from the disc using (3) from M. Jura 2003
    :param freq: The frequency of the light
    :param t_in: The inner temperature of the disc
    :param t_out: The outer temperature of the disc
    :param inc: The inclination of the disc
    :param wd: a WhiteDwarf object
    :return: The flux from the ring in mJy
    """
    def integrand(x):
        return x**(5/3)/np.expm1(x)

    x_in = H_P/(K_B*t_in)
    x_out = H_P/(K_B*t_out)

    A = 12*np.cbrt(np.pi)*(wd.radius/wd.distance)**2 * np.cos(inc) * ((2*K_B*wd.t_eff)/(3*H_P))**(8/3) * H_P/(C**2)

    f_ring = np.zeros(freq.size)

    for i in range(freq.size):
        f_ring[i] = A*np.cbrt(freq[i])*integrate.quad(integrand, np.abs(x_in*freq[i]), np.abs(x_out*freq[i]))[0]/1E-26

    return f_ring


def flux_int(freq, r_in, r_out, inc, wd):
    """

    :param freq: An array of frequencies in Hz.
    :param r_in: The inner radius of the disc in Solar radii.
    :param r_out: The outer radius of the disc in Solar radii.
    :param inc: The inclination of the disc in degrees.
    :param wd: A White Dwarf object.
    :return: An array of fluxes, in mJy
    """
    r = np.linspace(r_in, r_out, 4000)
    t_values = t_func(r, wd)
    r = r*R_SOLAR  # Conversion to cm after the t_values have been obtained

    def planck(nu, t):
        # The Planck function
        return 2 * H_P * nu**3/(C**2) * 1/np.expm1(H_P * nu/(K_B * t))

    fluxes = np.zeros(freq.size)
    const = 2 * np.pi * np.cos(inc)/wd.distance**2  # Flux integration constant.

    for i, nu in enumerate(freq):
        integrand = planck(nu, t_values) * r
        fluxes[i] = const * integrate.simps(integrand, r)

    return fluxes/1E-26

