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

def nh_ort1_make_stacks():
    
    """
    This program takes a directory full of individual NH KEM Hazard frames, stacks them, and stubracts
    a stack of background field. This reveals rings, etc. in the area.
    
    Written for NH MU69 ORT1, Jan-2018.
    
    """
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
    reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
    
    dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    hbt.figsize((12,12))
    
    # Load and stack the field images
    
    stack_field = image_stack(os.path.join(dir_data, reqid_field))
    flat_field  = stack_field.flatten(method = 'median')  # Median is much better than mean. Mean leaves CR's.
    
    for reqid in reqids_haz:
        stack_haz = image_stack(os.path.join(dir_data, reqid))
        flat_haz  = stack_haz.flatten(method = 'median')

        # Calculate the offset.
        # This is probably silly: we should get it from WCS instead
        
        out = ird.translation(flat_field, flat_haz)
        tvec = np.round(np.array(out['tvec'])).astype(int)
        
        # Make the plot
        
        plt.imshow(stretch(np.roll(np.roll(flat_haz,round(tvec[1]),1),round(tvec[0]),0) - flat_field))
        plt.title("{} - {}".format(reqid, "field"))
        plt.show()

if (__name__ == '__main__'):
    nh_ort1_make_stacks()
    