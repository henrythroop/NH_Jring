#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 22:36:55 2018

@author: throop
"""

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
import warnings
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
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)
import imreg_dft as ird                    # Image translation

from astropy.coordinates import SkyCoord

import os.path

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   get_radial_profile_circular import get_radial_profile_circular 

def nh_ort1_make_stacks():
    
    """
    This program takes a directory full of individual NH KEM Hazard frames, stacks them, and subtracts
    a stack of background field. This reveals rings, etc. in the area.
    
    Written for NH MU69 ORT1, Jan-2018.  
    """
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
    reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
    
    dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'

    zoom = 4
    
    # Set the edge padding large enough s.t. all output stacks will be the same size.
    # This value is easy to compute: loop over all stacks, and take max of stack.calc_padding()[0]
    
    padding = 61   
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
        
    # Set the RA/Dec of MU69. We could look this up from SPICE but it changes slowly, so just keep it fixed for now.
    
    radec_mu69 = (4.794979838984583, -0.3641418801015417)
    
    # Load and stack the field images
    
    stack_field = image_stack(os.path.join(dir_data, reqid_field))
    stack_field.align(method = 'wcs', center = radec_mu69)
    img_field  = stack_field.flatten(zoom=zoom, padding=padding)

    hbt.figsize((12,12))
    hbt.set_fontsize(20)
    
    for reqid in reqids_haz:
        stack_haz = image_stack(os.path.join(dir_data, reqid))
        stack_haz.align(method = 'wcs', center = radec_mu69)
        img_haz  = stack_haz.flatten(zoom=zoom, padding=padding)

        # Make the plot
        
        diff = img_haz - img_field
        diff_trim = hbt.trim_image(diff)
        plt.imshow(stretch(diff_trim))
        plt.title(f"{reqid} - field, zoom = {zoom}")
        plt.show()
        
        # Save the stacked image as a FITS file
        
        file_out = os.path.join(dir_data, reqid, "stack_n{}_z{}.fits".format(stack_haz.size[0], zoom))
        hdu = fits.PrimaryHDU(stretch(diff_trim))
        hdu.writeto(file_out, overwrite=True)
        print(f'Wrote: {file_out}')        

        # Make a radial profile
        
        pos =  np.array(np.shape(diff))/2
        (radius, profile) = get_radial_profile_circular(diff, pos=pos, width=1)
    
        hbt.figsize((10,8))
        hbt.set_fontsize(20)
        plt.plot(radius, profile)
        plt.xlim((0, 50*zoom))
        plt.ylim((-1,np.amax(profile)))
        plt.xlabel('Radius [pixels]')
        plt.title(f'Ring Radial Profile, {reqid}, zoom={zoom}')
        plt.ylabel('Median DN')
        plt.show()
        
if (__name__ == '__main__'):
    nh_ort1_make_stacks()
    