#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:39:08 2018

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
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astroquery.vo_conesearch import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
from tkinter import ttk
from tkinter import messagebox
tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

from   astropy.stats import sigma_clip

# HBT imports

import hbt

plt.set_cmap('plasma')

rj_km = 71492

# Now list all of the Gossamer observations. I have determined these groupings manually, based on timings.
# In general it looks like usually these are stacks of four frames at each pointing.

index_group = 6

index_image_list = [
#                    hbt.frange(46,48), # 1X1. Main ring. Maybe gossamer too but not interested.
#                    hbt.frange(49,53), # 1X1. Main ring. Maybe gossamer too but not interested.
#                    np.array([54]),    # 1X1. Main ring. Maybe gossamer too but not interested.
                    hbt.frange(59,62),  # Right ansa. I think Gossamer is there but hidden -- need to optimize.
                    hbt.frange(63,66),  # Left ansa. Furthest out, no ring. Use for subtraction.
                    hbt.frange(67,70),  # Left ansa. Far out, no ring. Closer than above.
                    hbt.frange(71,74),  # Left ansa. No ring. Closer than above.
                    hbt.frange(75,78),  # Left ansa. Gossamer. Closer than above.
#                    hbt.frange(79,88),  # Left ansa. Main ring. Prob gossamer too but not interested. 1X1. Svrl ptgs.
#                    hbt.frange(89,94),  # Io! Closest of the sequence, and I will ignore. 1X1.
                    hbt.frange(95,98),  # Right ansa. Gossamer limb.
                    hbt.frange(99,102), # Left ansa. Similar geometry as 75-78.
                    hbt.frange(112,115),# Right ansa. Gossamer limb. Similar geometry as 95-98 (offset pointing)
                    hbt.frange(116,119),# Left ansa. Similar geometry as 75-78. Sat visible??
                    hbt.frange(120,123),  # Left ansa. Closer than above.
                    hbt.frange(124,127),  # Left ansa. Closer than above.
                    hbt.frange(128,131),  # Left ansa + main ring ansa. Closer than above.
                    ]

#index_image = hbt.frange(71,74)
#index_image = hbt.frange(75,78)
#
#index_image = hbt.frange(95,99)
#
#index_image = hbt.frange(112,115)
#index_image = hbt.frange(116,119)
#index_image = hbt.frange(120,123)
#index_image = hbt.frange(124,127)
#index_image = hbt.frange(128,131)

#index_image = hbt.frange(100,102)
#index_image = hbt.frange(49,54)

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.


filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters 

dir_out    = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
dir_backplanes = '/Users/throop/data/NH_Jring/out/'

lun = open(dir_out + filename_save, 'rb')
t_all = pickle.load(lun)
lun.close()

abcorr = 'LT+S'

groups = astropy.table.unique(t_all, keys=(['Desc']))['Desc']

groupmask = t_all['Desc'] == groups[index_group]

t_group = t_all[groupmask]  # 

# =============================================================================
# Make a plot showing the center position of a lot of these frames.
# =============================================================================

file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel

sp.furnsh(file_tm)

hbt.figsize((10,5))

index_group = 6
index_images = np.array(
              list(hbt.frange(10,15)) +
              list(hbt.frange(46,54)) +
              list(hbt.frange(59,135)) )

#index_images = hbt.frange(49,54)

#index_images = np.array([46,47,48])

for index_images in index_image_list:  # Loop over *all* the images

    ra_arr      = []
    dec_arr     = []
    dist_rj_arr = []
    et_arr      = []
    dx_arr      = []
    dy_arr      = []
    dz_arr      = []
    corner_arr  = []

    arr_sum     = None

    plt.subplot(1,3,1)
    
#    fig, ax = plt.subplots()

    for index_image in index_images:   # Loop over the images in this obsevation
        
        t_i = t_group[index_image]  # Grab this, read-only, since we use it a lot.
        
        arr = hbt.read_lorri(t_i['Filename'])
        arr = hbt.lorri_destripe(arr)

        if arr_sum is None:
            arr_sum = arr
        else:    
            arr_sum += arr

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = WCS(t_i['Filename'])
    
        # Get the RA / Dec for the four corner points
        
        if t_i['Format'] == '1X1':
            radec_corners = w.wcs_pix2world([[0,0],[0,1024],[1023,1023],[1023,0],[0,0]],0) * hbt.d2r  # Radians
            radec_center  = w.wcs_pix2world([[511,511]],0) * hbt.d2r                       # Radians
        if t_i['Format'] == '4X4':
            radec_corners = w.wcs_pix2world([[0,0],[0,256],[256,256],[256,0],[0,0]],0) * hbt.d2r  # Radians
            radec_center  = w.wcs_pix2world([[128,128]],0) * hbt.d2r                       # Radians
        
    #    ra  = radec[:,0] 
    #    dec = radec[:,1]
    
        # Get RA / Dec for Jupiter
        
        (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'New Horizons')
        (_, ra_jup, dec_jup) = sp.recrad(vec[0:3])  # Radians
        vec_nh_jup = vec[0:3]
         
        plt.plot((radec_corners[:,0]-ra_jup)*hbt.r2d, (radec_corners[:,1]-dec_jup)*hbt.r2d,
                 label = f'{index_group}/{index_image}') # Plot all the four points for this box
        
        et_arr.append(t_i['ET'])
        
        is_limb_left = (radec_corners[0,0]-ra_jup)>0

    # Finished this loop over image set
    
    corner_arr = np.array(corner_arr) / rj_km  # This is now an array (n_pts x 5, 3)
    
    # Calculate Jupiter's radius, in degrees
    
    ang_jup_deg = (rj_km / sp.vnorm(vec_nh_jup)) * hbt.r2d
    
    circle1=plt.Circle((0,0),ang_jup_deg,   color='red') # Plot Jupiter itself
    circle2=plt.Circle((0,0),ang_jup_deg*2, facecolor='white', edgecolor='red') # Plot Jupiter itself
    circle3=plt.Circle((0,0),ang_jup_deg*3, facecolor='white', edgecolor='red') # Plot Jupiter itself
    
    plt.gca().add_artist(circle3)
    plt.gca().add_artist(circle2)
    plt.gca().add_artist(circle1)
    plt.gca().set_aspect('equal')
    plt.xlabel('Delta RA from Jupiter  [deg]')    
    plt.ylabel('Delta Dec from Jupiter [deg]')
    plt.xlim((4,-4))
    plt.ylim((-2,2))
    plt.legend(loc='upper right')
    plt.title(f'Gossamer images, {index_group}/{np.amin(index_images)}-{np.amax(index_images)}')

    # Print a message listing the range of R_J in this plot. 
    # **Or, label the x axis in RJ, on the plots themselves.**
    
    # 
    
    # Get the distance from Jupiter center to image center
    
    vec_nh_center = sp.radrec(1, radec_center[0][0], radec_center[0][1])  # Vector from NH, to center of LORRI frame
    ang_jup_center = sp.vsep(vec_nh_jup, vec_nh_center)                   # Ang Sep btwn Jup and LORRI, radians
    dist_jup_center_rj = ang_jup_center * sp.vnorm(vec_nh_jup) / rj_km    # Convert from radians into RJ
    width_lorri = 0.3*hbt.d2r                                             # LORRI full width, radians

    dist_jup_center_rj_range = np.array([ang_jup_center+width_lorri/2, ang_jup_center-width_lorri/2]) * \
        sp.vnorm(vec_nh_jup) / rj_km                                      # Finally, we have min and max dist
                                                                          # in the central row of array, in rj
    # NB: If we are on the RHS limb, then we need to reverse this.

    if not(is_limb_left):
        dist_jup_center_rj_range = dist_jup_center_rj_range[::-1]
    
    extent = [dist_jup_center_rj_range[0], dist_jup_center_rj_range[1],\
              dist_jup_center_rj_range[0], dist_jup_center_rj_range[1]]
    
    
         
#    plt.show() 
    
    # Plot the image, sfit subtracted
    
    plt.subplot(1,3,2)
    degree=5
    arr_sum /= len(index_images)  # We have just summed, so now divide
    arr_sum_s = hbt.remove_sfit(arr_sum,degree=degree)  # _s = sfit subtracted
    plt.imshow(stretch(arr_sum_s), origin='lower', extent=extent)
#    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.title(f'{index_group}/{np.amin(index_images)}-{np.amax(index_images)} - p{degree}')

    # Plot the image, with better subtraction
    
    plt.subplot(1,3,3)
    
    method = 'String'
    vars = '6/63-70 p5'
#    index_image = 124
#    index_group = 6
    arr_process = hbt.nh_jring_process_image(arr_sum, method, vars, 
                                     index_group=index_group, index_image=index_image, mask_sfit=None)


    plt.imshow(stretch(arr_process), origin='lower', extent=extent)
#    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.title(f'{index_group}/{np.amin(index_images)}-{np.amax(index_images)} - {vars}')

    # Now show all plots

    plt.tight_layout()
    plt.show()
    
    print(f'Is left: {is_limb_left}')
    print('---')

#%%%

# And hang on. We've already made backplanes for everything. Can I use those here? 
# I don't know. These images are obviously very edge-on.

#        file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')
#
#        # Save the shortname associated with the current backplane. 
#        # That lets us verify if the backplane for current image is indeed loaded.
#
#        file_backplane_shortname = t_group['Shortname'][index_image]
#				
#        if (os.path.isfile(file_backplane)):
#            lun = open(file_backplane, 'rb')
#            planes = pickle.load(lun)
#            lun.close()

# I could do RA/Dec relative to RA/Dec of Jupiter.
    