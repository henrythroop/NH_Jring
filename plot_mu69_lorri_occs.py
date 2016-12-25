#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 22:45:55 2016

@author: throop
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
import astropy.table as table
from astropy.coordinates import SkyCoord
import math
import hbt
import spiceypy as sp
    
file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69.tm'

# Read HD into a table. Its positions are 1950.

t = read_hd_all()

# Initialize 

sp.furnsh(file_tm)

utc_ca = '2019 1 jan 03:09:00'
et_ca  = sp.utc2et(utc_ca)

hour = 3600
day  = 24 * hour

num_dt = 1000

et = hbt.frange(et_ca + 0.1*hour, et_ca+1*day, num_dt)

ra_kbo = []
dec_kbo = []

for et_i in et:

    (st,lt) = sp.spkezr('MU69', et_i, 'J2000', 'LT+S', 'New Horizons')
    (range, ra_i, dec_i) = sp.recrad(st[0:3])
    
    ra_kbo.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo.append(dec_i)

ra_kbo = np.array(ra_kbo)
dec_kbo = np.array(dec_kbo)

plt.plot(ra_kbo*hbt.r2d, dec_kbo*hbt.r2d)
xlim = hbt.mm(ra_kbo*hbt.r2d)
ylim = hbt.mm(dec_kbo*hbt.r2d)

plt.plot(t['RA']*hbt.r2d, t['Dec']*hbt.r2d, linestyle='none', marker='.')
plt.xlim(xlim)
plt.ylim(ylim)
plt.show()

#plt.set_xlim(xlim)