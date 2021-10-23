# disc-models

This repository contains code I have written for the modelling of accretion discs. There are two main types of disc which are discussed below.
For each disc, the key feature is the temperature function. Once this has been obtained, we can model the spectrum of the disc using the Planck Function.

## Star.py

This class defines a stellar object. The main parameters are the object's radius,
effective temperature, and distance from Earth. When used with the disc modelling programs, the default units are solar masses, solar radii and kelvin
eg: if mystar.radius = 0.5 means 0.5 solar radii. 

In the equations below, properties of the star are represented with an asterisk in subscript. There are also functions to convert these to cgs units.

## flatmodel.py

This is an optically thick, geometrically thin disc described in M.Jura's 2003 Paper (https://arxiv.org/abs/astro-ph/0301411). 

t_func is the temperature function, equation (3) in the paper:

![equation](https://latex.codecogs.com/png.latex?%5Cbg_white%20T_%7Bring%7D%20%5Capprox%20%5Cleft%28%5Cfrac%7B3%7D%7B2%5Cpi%7D%20%5Cright%20%29%5E%7B1/4%7D%20%5Cleft%28%5Cfrac%7BR_*%7D%7BR%7D%5Cright%29%5E%7B3/4%7DT_*)

There are two functions to determine the flux from the disc. Both give return the same results, but one is more effecient.

flux_int creates an array of T values and R values, and uses these to integrate over the Planck function as shown in equation (2):

![equation](https://latex.codecogs.com/png.latex?%5Cbg_white%20F_%7Bring%7D%20%3D%20%5Cfrac%7B2%5Cpi%20cos%28i%29%7D%7BD_*%5E2%7D%5Cint_R_%7Bin%7D%5E%7BR_%7Bout%7D%7D%20B_%7B%5Cnu%7D%28T_%7Bring%7D%29RdR)

flux is the more effecient function, and uses a substitution to then integrate the function directly, shown in equation (4):

![equation](https://latex.codecogs.com/png.latex?%5Cbg_white%20F_%7Bring%7D%20%5Capprox%2012%5Cpi%5E%7B1/3%7D%5Cleft%28%5Cfrac%7BR_*%5E2%20cos%28i%29%7D%7BD_*%5E2%7D%20%5Cright%20%29%5Cleft%28%5Cfrac%7B2%20k_b%20T%7D%7B3h%5Cnu%7D%20%5Cright%29%5E%7B8/3%7D%5Cleft%28%5Cfrac%7Bh%5Cnu%5E3%7D%7Bc%5E2%7D%20%5Cright%20%29%5Cint_%7Bx_%7Bin%7D%7D%5E%7Bx_%7Bout%7D%7D%20%5Cfrac%7Bx%5E%7B5/3%7D%7D%7Be%5Ex%20-1%7Ddx)

where x = hv/kT

Both functions output the flux in terms of frequency, with units of mJy.

## flaredmodel.py

This models the optically thick but flared disc described in the E.Chiang et. al 2000 paper https://arxiv.org/abs/astro-ph/0009428.

t_eq uses Brent's root finder to solve equation A1 in appendix A of the paper:

![equation](https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Csin%5Cleft%5B%5Carctan%5Cleft%28%5Cgamma%5Cfrac%7BH%7D%7Bh%7D%5Csqrt%5Cfrac%7BT_i%7D%7BT_c%7D%5Csqrt%5Cfrac%7Ba%7D%7BR_*%7D%5Cright%29%20-%20%5Carctan%5Cleft%28%5Cfrac%7BH%7D%7Bh%7D%5Csqrt%5Cfrac%7BT_i%7D%7BT_c%7D%5Csqrt%5Cfrac%7Ba%7D%7BR_*%7D%20%5Cright%20%29%20&plus;%5Carcsin%5Cleft%28%5Cfrac%7B4%7D%7B3%5Cpi%7D%20%5Cfrac%7BR_*%7D%7Ba%7D%5Cright%20%29%20%5Cright%5D%20%3D%20%5Cfrac%7B2%7D%7B%5Cphi%7D%5Cleft%28%5Cfrac%7BT_i%7D%7BT_*%7D%20%5Cright%20%29%5E4%5Cleft%28%5Cfrac%7Ba%7D%7BR_*%7D%20%5Cright%20%29%5E2)

Presently, the opacity has not been included in the model.

get_t uses the method described in the paper to solve for T - it creates a logarithmic grid of r values, and solves for T at two consecutive values of r. It then updates the value of gamma.

Similarly to flux_int in the flatdisc.py, flux_int here uses the T and r values to integrate over the Planck function and determine the flux.
