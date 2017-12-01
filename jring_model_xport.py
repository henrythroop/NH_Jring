#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 19:52:42 2017

@author: throop
"""


import numpy as np
import astropy
import astropy.constants as c
import astropy.units     as u
import pytwobodyorbit


# mmu = gravitational parameter of centrla body. 1.3e20 for sun.
#     = G * m_central.
#       Units = m^3/s^2 (ie, mks)


r_jupiter = 71492*u.km
m_jupiter = 1.898e27*u.kg

mmu_jupiter = c.G * m_jupiter
a_adrastea = 129000*u.km
v_adrastea = np.sqrt(c.G * m_jupiter / a_adrastea).to('km/s')

mmu_jupiter = c.G * m_jupiter

dv = 10*u.m/u.s  # extra velocity, in mks

orbit = pytwobodyorbit.TwoBodyOrbit('Body', mmu=mmu_jupiter.value) # Convert to mks
pos = np.array( [a_adrastea.to('m').value, 0, 0] ) # Convert to mks
vel = np.array( [0, (v_adrastea + dv).to('m/s').value, 0] ) # Convert to mks

orbit.setOrbCart(0, pos, vel) # Triggers warnings, but ignore
kepl = orbit.elmKepl() 

xs, ys, zs, times = orbit.points(10000)
plt.plot(xs, ys, linestyle='none', marker='.')
plt.gca().set_aspect('equal')
plt.show()

# Histogram of x position (which is I guess projected profile)

(n, bins) = np.histogram(xs, bins=200, range = (0, np.amax(xs)))
plt.plot(bins[0:-1], n)
plt.xlabel('Projected X position')
plt.ylabel('Density')
plt.show()

# Make a histogram of radius. This *is* the radial profile. And it is 
# what we are after here. This is the goal.

rs = np.sqrt(xs**2 + ys**2 + zs**2)
(n, bins) = np.histogram(rs, bins=200, range = (np.min(rs), np.amax(rs)))
plt.plot((bins[0:-1]*u.meter / (1000*u.km)).decompose(), n)
plt.xlabel('Projected Radius [10^3 km]')
plt.ylabel('Density')
plt.title("J-ring gap, dv = {}".format(dv))
plt.show()
