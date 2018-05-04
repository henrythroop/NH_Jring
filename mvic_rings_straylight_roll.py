#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:15:29 2017

Program to do Q&D stray light analysis of MVIC Framing observations at a few different rotation angles.

This is intended for planning of the KEM MU69 flyby. We want to just add a few datapoints from the Pluto flyby
to the stray light vs. roll angle curve.

@author: throop
"""


# General python imports

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
from   subprocess import call
import warnings
import pdb
import os.path
import os
import subprocess

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling


import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

from astropy.convolution import Box1DKernel, Gaussian1DKernel, convolve
from   pymiecoated import Mie

from   scipy.stats import linregress

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure

# HBT imports

import hbt

# Local NH rings imports

#from  nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist

#from nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
#from nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image

dir_data = '/Users/throop/Data/NH_MVIC_Ring/'

file_tm = 'kernels/../gv_kernels_new_horizons.txt'  # SPICE metakernel

# Get a list of all the raw MVIC frames. We want to exclude Tod's mosaics, which are in form 'mvic...mos.fits'
# Also, we wants .fit file, not .fits. I'm not sure what the latter is.

files = glob.glob(dir_data + 'mpf*.fit')

# Loop over the files and get the sequence name

phase     = []
mean_dn   = []
median_dn = []
stdev_dn  = []
max_dn    = []
sequence  = []
et        = []
utc       = [] 
delta_et  = [] # Time from Pluto flyby, in seconds
name_sap  = []
exptime   = []
ang_phase = []
ang_roll  = []

# Start up SPICE

sp.furnsh(file_tm)

plt.set_cmap('Greys_r')

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

utc_pl = '2015 Jul 14 05:00:00'
et_pl = sp.utc2et(utc_pl)

for file in files:
    
    file_short = file.split('/')[-1].replace('.fit', '')
    hdulist = fits.open(file)
    shape = np.shape(hdulist['PRIMARY'].data)
    nframes = np.shape(hdulist['PRIMARY'])[0]
    
    # Grab the header. There is just one header total, and it has info on all the exposures.
    # e.g., DELTAT02 = mid-time for exposure #2.
    
    header = hdulist['PRIMARY'].header
    
    # Extract various quantities from the header

    et_i       = header['SPCSCET']
    delta_et_i = et_i - et_pl
    utc_i      = sp.et2utc(et_i, 'C', 0)
    name_sap_i = header['SAPNAME']
    exptime_i  = header['EXPTIME']
    
    # Make sure the SAP name has 'RING' or 'ring' in it somewhere
    
    if 'RING' in name_sap_i.upper():
        et_i = header['SPCSCET']
        
        ang_roll_i = 0
        if name_sap_i == "O_RingDep_A_1":
            ang_roll_i = 180
        else:
            ang_roll_i = 45
        
        delta_et.append(delta_et_i)
        et.append(et_i)
        utc.append(utc_i)
        exptime.append(exptime_i)
        name_sap.append(name_sap_i)
        ang_roll.append(ang_roll_i)
        
        print ("File {}, {}, K+ {:3.1f} d, {}, {}".format(file, utc_i, delta_et_i/(3600 * 24), name_sap_i, shape))
        
        # Measure some quantities from the image itself
        
    #    for nframe in range(nframes):
    
        nframe = 1 # Each FITS file can have a stack of many frames. We just do stats on one. Which one?
        
        data = hdulist['PRIMARY'].data[nframe,:,:]
        median_dn.append(np.median(data))
        mean_dn.append(np.mean(data))
        stdev_dn.append(np.std(data))
        
        et_i = et[-1]  # Get the value we just pushed there    
        (st_nh_pl, lt)  = sp.spkezr('Pluto', et_i, 'J2000', 'LT+S', 'New Horizons')
        (st_sun_nh, lt) = sp.spkezr('New Horizons', et_i, 'J2000', 'LT+S', 'Sun')
        vec_nh_pl  = st_nh_pl[0:3]
        vec_sun_nh = st_sun_nh[0:3]
        phase_i = sp.vsep(vec_sun_nh, vec_nh_pl) * hbt.r2d
        phase_i = np.around(phase_i,decimals=1)   # Round the angle (55.61987 -> 55.6)
        ang_phase.append(phase_i)

        # View the image on the screen
        
        ratio_aspect = 3
        
        hbt.set_fontsize(18)
        hbt.figsize((15,9))
        plt.subplot(3,1,1)
        plt.imshow(stretch(data), aspect=ratio_aspect)
        plt.ylabel('Y [pix]')
#        plt.xlabel('X [pix]')
        plt.title("{}, {}, Raw".format(name_sap_i, file_short))
        plt.show()

        plt.subplot(3,1,2)
        plt.imshow(stretch(hbt.remove_sfit(data)), aspect=ratio_aspect)
        plt.ylabel('Y [pix]')
#        plt.xlabel('X [pix]')
        plt.title('Raw - Polynomial')
        plt.show()

        plt.subplot(3,1,2)        
        cols = np.median(data, axis=0)/exptime_i
        plt.plot(cols)
        plt.title('Raw, Column Median')
        plt.ylabel('DN/sec')
        plt.xlabel('X [pix]')
        plt.show()
        
# Ideally we'd put all these values into a table, or at least some NP arrays. Haven't done that yet, though.

# Make some plots

ang_roll  = np.array(ang_roll)
exptime   = np.array(exptime)
median_dn = np.array(median_dn)
mean_dn   = np.array(mean_dn)
stdev_dn  = np.array(stdev_dn)
        
plt.plot(ang_phase, median_dn, marker = '.', linestyle = 'none', label = 'DN Median')
plt.plot(ang_phase, mean_dn, marker = '.', linestyle = 'none', label = 'DN Mean')
plt.plot(ang_phase, stdev_dn, marker = '.', linestyle = 'none', label = 'DN Stdev')

plt.ylabel('DN')
plt.xlabel('Phase angle')
plt.show()

#plt.plot(et)
plt.show()

ms = 16
alpha = 0.4

hbt.set_fontsize(20)
plt.plot(ang_roll, median_dn/exptime, marker = '.', linestyle = 'none', label = 'DN/sec Median', ms=ms, alpha=alpha)
plt.plot(ang_roll, mean_dn/exptime, marker = '.', linestyle = 'none', label = 'DN/sec Mean', ms=ms, alpha=alpha)
#plt.plot(ang_roll, stdev_dn/exptime, marker = '.', linestyle = 'none', label = 'DN Stdev')
plt.ylabel('DN/sec')
plt.legend()
plt.xlabel('Roll Angle from SC +Y [deg]')
plt.show()

    # Compute the roll angle and phase angle
   

#def get_angle_sun():
#
##; Calculates the projected angle to the sun, from view of the observer.
##; This is useful for stray light calculations, so we can see direction to the sun,
##; even if it is off the screen.
##; 
##; This is lightly modified from GV_GET_ANGLE_PLOT_ROTATE.
##; Both routines are quite a bit more complex than one would think necessary.
##
##; HBT 1-Nov-2017
#
#; Returns position angle of Sun, CCW, from NCP, in radians.
#
#common units
#common gv_common
#
#  index_dt = 0
#
#  vec_x_sc = np.array([1., 0., 0.])		# S/C +X vector
#  vec_y_sc = np.array([0., 1., 0.])		# S/C +Y vector, good ref for NH
#  vec_z_sc = np.array([0., 0d, 1.])		# S/C +Z vector
#
#  vec_north_j2k = np.array([0., 0., 1.])		# NCP or NEP -- should work for either one.
#
#  name_body_angle = 'Sun'		; Which body do we do it relative to. Usually Sun.
#
#  index_body	= (where(strupcase(name_bodies) eq strupcase(name_body_angle), found))[0]
#  vec_obs_body   = (eph[index_body, index_dt].state)[0:2]
#
## Get the s/c boresight vector
#
#    'NEW HORIZONS' : vec_boresight_sc	= -vec_x_sc		# NH: Boresight is x
# 
## Define a plane which is perpendicular to the S/C boresight. We will project sun vector into this plane.
## This plane is defined in S/C frame. That makes the math particularly easy, because while we will end up with 
## 3D vectors, they will actually be just 2D when projected.
#
#      CSPICE_NVC2PL, vec_boresight_sc, 0d, plane_sc		; "Normal vector and constant to plane"
#
## Calculate the matrix to convert from J2K to s/c frame
## NB: This assumes that the rotation is defined by mx_inst_frame. That is true in case of SPICE lookup, but not if the 
## user has pointed the
## FOV manually. In that case, we can look it up this way, as per gv_plot_fovs.pro:
#    
#      if (name_target ne 'C-KERNEL') then begin
#
#        mx_inst_frame	= GV_GET_FOV_ROTATE_MX(ra_lon_fov_center, dec_lat_fov_center, $		
#                               fov_rotate, vec_boresight_sc, $
#			       vec_axis_scan=vec_axis_scan, $
#			       vec_ext=vec_ext)	
#      end
#
#      CSPICE_INVERT, mx_inst_frame, mx_frame_inst
#
## Rotate North pole into s/c frame, based on the rotation matrix from SPICE lookup
#
#      vec_north_sc = GV_MXVG(mx_frame_inst, vec_north_j2k)
#
## Rotate observer-to-Sun vector into s/c frame, based on the rotation matrix from SPICE lookup
#
#      vec_body_sc = GV_MXVG(mx_frame_inst, vec_obs_body)
#
## Now project sun vector into this plane
## Project vector to body into plane perp to s/c boresight. Trivial.
#
#      CSPICE_VPRJP, vec_body_sc, plane_sc, vec_body_sc_prj	
#
## And project NCP into this plane
#
#      CSPICE_VPRJP, vec_north_sc, plane_sc, vec_north_sc_prj	; Project NCP into plane perp to s/c boresight.
#
# Now take angle from NCP, to S/C X.
# We can't use CSPICE_VSEP, since that always returns value in 0..pi, *and* returns same angle whether CCW or CW rotation. 
# I think SPICE should have a function like VSEP_SIGNED. So we have to do it manually, which took a lot of 
# effort to figure out.
# See http://d-rob.org/blog/2011/06/angle-vector-to-another/ . I use his final 
# eq (projected into 2D), not his penultimate one (which is wrong).
#
#      v1 = vec_north_sc_prj
#      v2 = vec_body_sc_prj
#
#      sep      = atan(v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1])  ; Only for 2D, not 3D. 
#      # http://d-rob.org/blog/2011/06/angle-vector-to-another/
#
#      if ((v1[0] eq 0) and (v2[0] eq 0)) then $
#        sep      = atan(v1[1]*v2[2] - v1[2]*v2[1], v1[1]*v2[1] + v1[2]*v2[2])  ; Same as above
#
#      ; Convert it into a clockwise angle
#
#      sep = -sep
#
#      ; Make it a positive angle
#
#      if (sep < 0) then sep = sep + 2d * pi
#
#  return, sep  ; Return angle in radians, CW
#
#end

    
    