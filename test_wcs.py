#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 14:20:36 2018

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

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular

from   get_radial_profile_circular import get_radial_profile_circular
# HBT imports

import hbt
from image_stack import image_stack
from plot_img_wcs import plot_img_wcs
from wcs_translate_pix import wcs_translate_pix 
    
# =============================================================================
# Start of code
# =============================================================================
   
stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

plt.set_cmap('Greys_r')

zoom = 1      # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
              # 1 is faster; 4 is slower but better.

#    name_ort = 'ORT1'
#    name_ort = 'ORT2_OPNAV'
              
name_ort = 'ORT2'
initials_user = 'HBT'
dir_data = '/Users/throop/Data'

dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
dir_out       = os.path.join(dir_data, name_ort)
reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
reqids_haz  = ['K1LR_HAZ02', 'K1LR_HAZ03']
reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'

# Start up SPICE if needed

if (sp.ktotal('ALL') == 0):
    sp.furnsh('kernels_kem_prime.tm')

hbt.figsize((12,12))
do_force = False   # If set, reload the stacks from individual frames, rather than restoring from a pickle file.

# Load and stack the field images, one stack.

stack_field = image_stack(os.path.join(dir_images, reqid_field),   do_force=do_force, do_save=do_force)

# Plot some frames to check the WCS. Everything looks OK.

plot_img_wcs(  (stack_field.t[0]['data']), stack_field.t[0]['wcs'])
plot_img_wcs(  (stack_field.t[1]['data']), stack_field.t[1]['wcs'])
plot_img_wcs(  (stack_field.t[2]['data']), stack_field.t[2]['wcs'])
plot_img_wcs(  (stack_field.t[2]['data']), stack_field.t[2]['wcs'])
plot_img_wcs(  (stack_field.t[2]['data']), stack_field.t[2]['wcs'])
        
# Try doing a wcs translation on one image. Shift by fixed # of pixels.

dx_pix = -40  # -10, -40 is the actual shift for field[0]
dy_pix = -50

pad    = 50
pad_xy = [[pad,pad], [pad,pad]]

# Get field WCS

wcs = stack_field.t[0]['wcs']
img = stack_field.t[0]['data']
plot_img_wcs(img, wcs, title = 'Field 0')  # Plots the image, with a point indicated

# Copy field WCS into new WCS2.

wcs2 = wcs.deepcopy()

# Apply translation to offset the WCS

wcs_translate_pix(wcs2, -dy_pix-pad, -dx_pix-pad)

# Apply roll to offset the image

img2 = np.roll(np.roll(np.pad(img, pad_xy, mode='constant'), dx_pix, axis=0), dy_pix, axis=1)

# Now plot and see that they match. [A: Yes, it works, for many different combinations of dx dy.]

plot_img_wcs(img2, wcs2, title = f'Rolled dx = {dx_pix}, dy= {dy_pix}')  # Confirmed: this works. So, the basics of 
                                                                        # my image translation work.

# Look up the position of MU69

radec_mu69 = (4.794984626030024, -0.364179487109049) # Keep in radians

# Align the field frames on a specific RA/Dec. What this does is set the offset value in the stack
# Output is visible in stack_field.t['shift_x_pix'][0]

stack_field.align(method = 'wcs', center = (radec_mu69))

# Align Hazard frames on specifc RA/Dec

#for reqid_i in reqids_haz:
#    stack_haz[reqid_i].align(method  = 'wcs', center = (radec_mu69))  # In each individual stack, align all images

# Calc the padding required. This can only be done after the images are loaded and aligned.

pad_field = stack_field.calc_padding()[0]
pad_haz = []
for reqid_i in reqids_haz:
    pad_haz.append(stack_haz[reqid_i].calc_padding()[0])
pad_haz.append(pad_field)
    
pad = max(pad_haz)

# Flatten the stacks into single output images

# Flatten the field stack
# This should return a correct WCS. But it does not! This is right here the problem.
# XX stack_field.t['shift_x_pix'][0] = -10
# XX stack_field.t['shift_y_pix'][0] = -40


(img_field, wcs_field) = stack_field.flatten(do_subpixel=False, method='median',zoom=zoom, padding=pad)

# Flatten the main hazard stacks
# When we do this, the output image is shifted around within the padding amount, so that *all* 
# ouput images have the same size. So, maybe flatten should return img *and* wcs?
# And the expectation would be that all of these individual WCS would match. That would be the point.

img_haz = {}
wcs_haz = {}
img_haz_diff = {}

for reqid_i in reqids_haz:
    (img_haz[reqid_i], wcs_haz[reqid_i])  =\
          stack_haz[reqid_i].flatten(do_subpixel=False,  method='median',zoom=zoom, padding=pad)
    img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field
    
    plot_img_wcs(img_haz[reqid_i], wcs_haz[reqid_i], title = reqid_i)

