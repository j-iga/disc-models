"""
A program containing functions useful for the modelling of flared accretion discs.
Uses the results of E.Chiang et al. 2001. https://arxiv.org/abs/astro-ph/0009428

Author: Jordan Iga, The University of Edinburgh
        Imperial College London UROP Student
"""

from scipy import integrate
from scipy import optimize
import numpy as np
from Star import Star

# UNIVERSAL CONSTANTS https://physics.nist.gov/cuu/Constants/
C = 29979245800  # Speed of light in a vacuum (cm/s)
H_P = 6.62607015E-27  # Planck constant in cm^2 g/s
K_B = 1.380649E-16  # Boltzmann constant in cm^2 g/s^2/K
G = 6.67430E-8  # Newton's constant of Gravitation in cm^3/s^2
M_SOLAR = 1.989E33  # Solar Mass in grams
PARSEC = 3.0857E18  # Parsec in cm
R_SOLAR = 6.957E10  # Solar radius in cm
AU = 1.49598073E13  # Astronomical Unit in cm


def t_eq(r, g, h, mu, phi, wd):
    """
    Equation (A1) from E.Chiang et al. 2001

    :param r: The distance from the WD.
    :param g: Gamma, (A2) from E.Chiang et al. 2001.
    :param h: Visible photospheric height/gas scale height.
    :param mu: The mean molecular density of the disc.
    :param phi: A constant that is approximately 0.5
    :param wd: A Star object.
    :return:
    """
    #t = r
    #r = r*R_SOLAR
    t_c = G*Star.mass_cgs(wd)*mu/(K_B*Star.radius_cgs(wd))

    def f_t(t):
        f = np.sin(np.arctan(g*h*np.sqrt(t/t_c)*np.sqrt(r/wd.radius)) - np.arctan(h*np.sqrt(t/t_c)*np.sqrt(r/wd.radius))
                   +np.arcsin(4/(3*np.pi) *wd.radius/r)) - (2/phi)*(t/wd.t_eff)**4 * (r/wd.radius)**2
        return f

    temp = optimize.brentq(f_t, 0, 3000)

    return temp


def get_t(r_in, r_out, gamma0, h_ratio, mu, phi, wd):
    """
    Implements the method of solving for T as described in the appendix of
    E.Chiang et al. 2001.

    :param r_in: The inner radius of the disc.
    :param r_out: The outer radius of the disc.
    :param gamma0: The initial guess for gamma.
    :param h_ratio: Visible photospheric height/gas scale height.
    :param mu: The mean molecular mass of the disc.
    :param phi: A constant that is approximately 0.5.
    :param wd: A Star object.
    :return: An array of temperatures for the disc, and an array of their respective distances.
    """
    N = 400 # Number of distance steps
    gamma_i = gamma0

    radii = np.logspace(np.log10(r_in), np.log10(r_out), N)
    t_values = np.zeros(N)

    for i in range(N):
        if i % 2 == 0:
            t_i1 = t_eq(radii[i], gamma_i, h_ratio, mu, phi, wd)
            t_i2 = t_eq(radii[i+1], gamma_i, h_ratio, mu, phi, wd)
            t_values[i] = t_i1
            t_values[i+1] = t_i2

            gamma_i = (3/2) + (1 / 2) * np.log(t_i2 / t_i1) / np.log(radii[i+1] / radii[i])

    return radii, t_values


def flux_int(r_in, r_out, freq, gamma, h_ratio, mu, phi, wd):

    r, t_values = get_t(r_in, r_out, gamma, h_ratio, mu, phi, wd)
    r = r*R_SOLAR

    def planck(nu, t):
        # The Planck function
        return 2 * H_P * nu**3/(C**2) * 1/np.expm1(H_P * nu/(K_B * t))

    fluxes = np.zeros(freq.size)

    const = 2 * np.pi / Star.distance_cgs(wd)**2

    for i, nu in enumerate(freq):
        integrand = planck(nu, t_values) * r
        fluxes[i] = const * integrate.simps(integrand, r)

    return fluxes/1E-26
