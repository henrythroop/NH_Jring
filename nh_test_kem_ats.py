#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

NH_TEST_KEM_ATS -- test some trajectory files

This is just a simple one-off program to do some tests on YPG's KEM flyby trajectories.

Concl of this problem: She used one NAIF ID for MU69, and the earlier trajectory I had was using a different one.

Created on Tue Jun 27 22:14:44 2017

@author: throop
"""

import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
import os.path

from   astropy.utils import data
from   astropy.wcs import WCS

import astropy
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np

import spiceypy as sp

# HBT imports

import hbt

sp.furnsh('kernels_kem_ats.tm')

utc_start = "2019 1 Jan 06:40:00"
utc_end   = "2019 1 Jan 07:20:00"

et_start = sp.utc2et(utc_start)
et_end   = sp.utc2et(utc_end)

et_mid = (et_start + et_end)/2

num_dt = 500

et = hbt.frange(et_start, et_end, num_dt)

name_target = 'MU69'
name_observer = 'New Horizons'

dist   = np.zeros(num_dt)
phase  = np.zeros(num_dt)

for i,et_i in enumerate(et):
    (state,_) = sp.spkezr(name_target, et_i, 'J2000', 'LT+S', name_observer)
    dist[i] = np.sqrt(np.sum(state[0:3]**2))
    
plt.plot(et - et_mid, dist)
plt.xlabel('Dist [km]')
plt.ylabel('Seconds from c/a')
plt.show()    
    
print("dist[0] = {}".format(dist))
