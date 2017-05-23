#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:11:33 2017

This is a one-off code to convert a table of METs to phase angle, for Tod Lauer
The datafile is read from a draft I extracted from 'v4' of the ring search paper.
The goal here is to add a 'phase angle' column to Table 1.

@author: throop
"""

from   astropy.table import Table
import numpy as np
import spiceypy as sp
import hbt

file_tm = 'kernels/../gv_kernels_new_horizons.txt'  # SPICE metakernel
file_met = 'list_met.txt'
file_met_phase = 'list_met_phase.txt'

sp.furnsh(file_tm)

table     = Table.read(file_met, format = 'ascii')      # Read Tod's table from disk
met       = table['MET'].data
phase_deg = np.zeros(np.shape(met))

phase_list = []

for met_i in met:
    utc_i = hbt.met2utc(met_i)
    et_i  = sp.utc2et(utc_i)
    
    (st_nh_pl, lt)  = sp.spkezr('Pluto', et_i, 'J2000', 'LT+S', 'New Horizons')
    (st_sun_nh, lt) = sp.spkezr('New Horizons', et_i, 'J2000', 'LT+S', 'Sun')
    
    vec_nh_pl  = st_nh_pl[0:3]
    vec_sun_nh = st_sun_nh[0:3]
    
    phase_i = sp.vsep(vec_sun_nh, vec_nh_pl) * hbt.r2d
    phase_i = np.around(phase_i,decimals=1)   # Round the angle (55.61987 -> 55.6)
    phase_list.append(phase_i)

table['Phase [deg]'] = np.array(phase_list)  # Save into an output table

table.write(file_met_phase, format = 'ascii', overwrite=True)
print("Wrote: " + file_met_phase)            # Write back to disk, in new file

