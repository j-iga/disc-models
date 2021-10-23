# disc-models

This repository contains code I have written for the modelling of accretion discs. There are two main types of disc which are discussed below.
For each disc, the key feature is the temperature function. Once this has been obtained, we can model the spectrum of the disc using the Planck Function.

## flatmodel.py

This is an optically thick, geometrically thin disc described in M.Jura's 2003 Paper (https://arxiv.org/abs/astro-ph/0301411). The temperature T of the disc at 
a radius R is given by 

![equation](https://latex.codecogs.com/svg.latex?T_%7Bring%7D%20%5Capprox%20%5Cleft%28%5Cfrac%7B3%7D%7B2%5Cpi%7D%20%5Cright%20%29%5E%7B1/4%7D%20%5Cleft%28%5Cfrac%7BR_*%7D%7BR%7D%5Cright%29%5E%7B3/4%7DT_*)
