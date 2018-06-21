#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 10:09:56 2018

@author: throop
"""

import glob
import math
import os.path
import os

import astropy
from   astropy.io import fits
import astropy.table
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library
import pickle # For load/save

import scipy

# HBT imports

import hbt

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular
from   plot_img_wcs import plot_img_wcs
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes

# =============================================================================
# 
# Idea here is to plot a ring surrounding MU69. This code will iterate over many 
# possibilities for RA/Dec for the pole vector, in order to find the right one.
# 
# - Read the FITS file from the superstack. A small extract of it.
# - Plot that to screen.
# - On top of that, plot a ring centered on MU69, as seen from NH.
# - That ring is defined by a plane and a point.
#     - Plane is normal to a specified RA/Dec
#     - Point is center of UT, in J2000 space
# - Drop 100 points onto this plane, and plot their location.  
# =============================================================================

def nh_ort_find_ring_pole():
    
    file_superstack = '/Users/throop/Data/ORT4/superstack_ORT4_z4_mean_wcs_sm_hbt.fits'
    
    file_tm = 'kernels_kem_prime.tm'
    sp.unload(file_tm)
    sp.furnsh(file_tm)
    
    f = fits.open(file_superstack)
    
    img = f[0].data
    
#    plt.imshow(stretch(img))
#    plt.show()
    
    wcs = WCS(file) 
    
    num_pts = 200
    
    ra_pole   = 275 * hbt.d2r
#    dec_pole  = -56 * hbt.d2r
    dec_pole  = 13 * hbt.d2r

    radius_ring = 9000  # Radius in km

    vec_pole_j2k = sp.radrec(1, ra_pole, dec_pole)
    
    et = float(f[0].header['SPCSCET'])
    utc = sp.et2utc(et, 'C', 0)
    
    # Get position from NH to UT
    
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')

    vec_nh_ut = st[0:3]
 
     # Get position from Sun, to NH
    
    (st, lt) = sp.spkezr('New Horizons', et, 'J2000', 'LT', 'Sun')
    
    vec_sun_nh = st[0:3]

    vec_sun_ut = vec_sun_nh + vec_nh_ut
    
    # Define a 'ring plane', based on a pole vector, and a point
    # This ring plane should be in J2K space -- that is, centered on Sun.
    
    plane_ring = sp.nvp2pl(vec_pole_j2k, vec_sun_ut)  # Pole position is variable. Point is UT in J2K.
 
    # Get the point and spanning vectors that define this plane
    
    # XXX for some reason, these values from pl2psv do not depend on value of vec_pol_j2k
    
    (pt_pl, vec1_pl, vec2_pl) = sp.pl2psv(plane_ring)

    # Now take a bunch of linear combinations of these spanning vectors
    
    # Plot UT's position on the plot
    
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')

    (_, ra, dec) = sp.recrad(vec_nh_ut)
        
    (x, y) = wcs.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
       
#    plt.plot(x, y, marker = 'o', ms = 10, alpha=0.3, color='purple')
    
    # Set an offset to the WCS values, in case UT is not in the right position (ie, not centered properly)
    
    dy = 0 # Large value moves up
    dx = 0

    # Draw the ring image
    
    plt.imshow(stretch(img), origin='lower')
    
    # Calculate and draw all of the ring points
        
    for i in range(num_pts):
        
        angle_azimuth = 2*math.pi * (i / num_pts)   # Put in range 0 .. 2 pi
        vec_i = vec1_pl * math.sin(angle_azimuth) + vec2_pl * math.cos(angle_azimuth)
        vec_i = vec_i * radius_ring
        
        # Now get the point in space, J2K
        
        pt_ring_i_j2k = vec_i + vec_sun_ut
        
        vec_sun_ring_i = pt_ring_i_j2k
        
        vec_nh_ring_i = vec_sun_ring_i- vec_sun_nh
        
        # 
        (_, ra_i, dec_i) = sp.recrad(vec_nh_ring_i)
        
        (x, y) = wcs.wcs_world2pix(ra_i*hbt.r2d, dec_i*hbt.r2d, 0)
        
        plt.plot(x+dx, y+dy, marker = 'o', ms = 1, color='red', alpha = 0.15)
        print(f'{i}, {ra_i*hbt.r2d}, {dec_i*hbt.r2d}, {x}, {y}')
    
    plt.title(f'ORT4 Superstack, Ring Pole = ({ra_pole*hbt.r2d},{dec_pole*hbt.r2d}) deg')    
    plt.show()
    
    return

if (__name__ == '__main__'):
    nh_ort_find_ring_pole()
    
    