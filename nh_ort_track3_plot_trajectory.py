#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 14:38:05 2018

@author: throop
"""

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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)
import imreg_dft as ird                    # Image translation
import struct

import re # Regexp
import pickle # For load/save

from   datetime import datetime

import scipy

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular

# HBT imports

import hbt

# =============================================================================
# Do a one-off test to make some plots to validate the XYZ orientation.
# This does not merge or do anything else. It just makes nicely scaled plots
# of the Track-3 output, to visualize the trajectories and disk orientations.
#    
# NB: Large beta → Small grains    
# =============================================================================
    
def nh_ort_track3_plot_trajectory():

    plt.set_cmap('plasma') 
    dir =  '/Users/throop/data/ORT2/hamilton/deliveries/'
    runs_full = [ dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta1.0e-04/subset07',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003_original/y3.0/beta1.0e-04/subset07',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta1.0e-03/subset02',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003_original/y3.0/beta1.0e-03/subset02',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta2.2e-02/subset02',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta2.2e-02/subset05']

    # For the order, we always want (0,0) to be lower left.
    
    origin = 'lower'                               # Does imshow() plot start w/ (0,0) at lower-left, or upper-left?

    # Define the axes. DK says that the array is delivered in order [x, y, z], which are same as Mark's coord frame.
    # That means that if I sum in the '0' direction, I will have a plot in Y and Z.
    
    axes           = ['X',     'Y',     'Z']    
    num_axis       = {'X' : 0, 'Y' : 1, 'Z' : 2}
    
    # Now, make a dictionary to show us what the axes will be of the output image after doing np.sum().
    # The dictionary here means: 
    #     If we sum along the X axis (0), and we plot the result, what will be on vertical axis of the imshow() plot.
    #     In this case, first remaining axis (Y) is on vertical, and second (Z) is on horizontal.
    #     That is how imshow() works.
    
    axes_vertical  = {'X':'Y', 'Y':'X', 'Z':'X'}   
    axes_horizontal= {'X':'Z', 'Y':'Z', 'Z':'Y'}
    
    hbt.figsize((18,10))
    hbt.fontsize(10)  # Default
    fontsize_axes = 15
    
    do_stretch_linear = False
    
    for run_full in runs_full:
        i = 1
        run  = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'           
        ring = nh_ort_track3_read(run)
        ring.print_info()
    
        for axis in axes:

            plt.subplot(1,3,i)

            # Plot the individual image. Start it with origin at lower-left corner.

            img = np.sum(ring.density, axis=num_axis[axis])  # Make the flattened image
            width_colorbar = 10
            img_stretch = stretch_hbt(img)
            if do_stretch_linear:
                img_stretch = astropy.visualization.PercentileInterval(98)(img)
            
            # Create the colorbar, and superimpose it on the image
            
            colorbar = hbt.frange(np.amin(img_stretch), np.amax(img_stretch), hbt.sizey(img_stretch))
            img_stretch[:,-1] = colorbar  # There is probably a better way to do this?
            img_stretch[:,-2] = colorbar
            img_stretch[:,-3] = colorbar
            img_stretch[:,-4] = colorbar
            img_stretch[:,-5] = colorbar
            
            # Display the image.
            # The image is displayed in exactly the same orientation as if I print it, with the exception
            # that the origin={lower | upper} keyword can flip it vertically.
            # When accessing the array elements, they are in order img[y, x] -- which is opposite IDL.
            
            plt.imshow(img_stretch, extent=extent, origin=origin)
            
            # Create the labels for the colorbar, and individually place them
            
            num_ticks_colorbar = 5 # Number of vertical value to put on our colorbar
            
            for j in range(num_ticks_colorbar):
                val = stretch_invert(hbt.frange(np.amin(img_stretch), np.amax(img_stretch), num_ticks_colorbar)[j])
                val = round(val)  # Convert to zero if it is very close
                xval = 0.65*np.max(extent)
                yrange = np.max(extent)-np.min(extent)
                yval = np.min(extent) + 0.02*yrange + (j/(num_ticks_colorbar-1) * 0.92*yrange)
                if not(do_stretch_linear):
                    plt.text(xval, yval, f'{val:.1e}', color = 'white') # Label the colorbar, if log stretch only.
            
            # Label the axes and the plot
            
            plt.title(f'Summed along {axis}', fontsize=fontsize_axes)
            plt.ylabel(axes_vertical[axis] + ' [km]', fontsize=fontsize_axes)
            plt.xlabel(axes_horizontal[axis] + ' [km]', fontsize=fontsize_axes)
            plt.tight_layout()
            
            i+=1
        plt.show()
        print('\n-----\n')
        

if (__name__ == '__main__'):
    nh_ort_track3_plot_trajectory()