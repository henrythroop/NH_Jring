#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:32:23 2017

@author: throop
"""

# General python imports

import pdb
import glob
import warnings
import os.path
import os

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
from   astropy.utils import data

import spiceypy as sp
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
from   scipy.stats import mode
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   photutils import DAOStarFinder
#import cv2

import re # Regexp
import pickle # For load/save

from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure
import warnings
from   importlib import reload
from   time import gmtime, strftime

#from  nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist

# HBT imports

import hbt

#==============================================================================
# NH_JRING_MASK_FROM_OBJECTLIST
# 
# Given the name of an objectlist, creates a boolean array mask with one value per pixel,
# indicating positions of stars, satellites, etc according to WCS, SPICE, and star catalogs.
#==============================================================================

def nh_jring_mask_from_objectlist(objectfile, do_plot = False):
    
    '''
    Input: an objectfile ('_opnav_objects.txt')
    Output: an array with a True/False mask. True = object there. False = no object.
    
    Boolean array mask with one value per pixel, indicating positions of stars, 
    satellites, etc according to WCS, SPICE, and star catalogs.
    '''

# Define the width of the 'PSF' for each type of object
    
    width_star     = 10
    width_sat      = 25

    width_star_4x4 = 3
    width_sat_4x4  = 8
    
# Set the directories and paths
    
#    objectfile = 'lor_0034604523_0x630_sci_1_opnav_objects.txt'
    
    dir_images = '/Users/throop/Data/NH_Jring/data/jupiter/level2/lor/all/'
    dir_out    = '/Users/throop/Data/NH_Jring/out/'
    
    objectfile_base = objectfile.split('/')[-1]
    path_objectfile = dir_out + objectfile_base

# Extract the original image file (the opnav'd FITS file)
    
    file_image = objectfile_base.replace('_objects', '').replace('.txt', '.fit')
    path_image = dir_images + file_image
    
    header = hbt.get_image_header(path_image)
    dx_pix = header['NAXIS1']
    dy_pix = header['NAXIS2']
    mode   = header['SFORMAT']
    
    t = Table.read(path_objectfile, format = 'csv')

# Create the mask functions

    if (mode == '4X4'):
        width_star = width_star_4x4
        width_sat  = width_sat_4x4
        
    mask_star = (hbt.dist_center(width_star, centered = True, invert=False) <= width_star/2.)    
    mask_sat  = (hbt.dist_center(width_sat , centered = True, invert=False) <= width_sat /2.)    

# Create the output array

    mask = np.zeros((dx_pix, dy_pix))

    xx, yy = np.mgrid[:dx_pix, :dy_pix]

# Loop and create an entry for each object
    
    for row in t:
        x = row['x_pix']
        y = row['y_pix']
        name = row['name']
#        print("{}: x={}, y={}".format(name, x, y))
        
        if (name == 'star'):
            width = width_star
        else:
            width = width_sat
            
        mask_i = np.sqrt( (xx-x)**2 + (yy-y)**2 ) < width
        
        mask = np.logical_or(mask, mask_i)

# Plot it if requested

    if (do_plot):
        plt.imshow(mask)
        plt.title("{}, N = {}".format(objectfile_base, np.size(t)))
        plt.show()

# Transpose it, which happens to be what it needs for LORRI images

    mask = np.transpose(mask)
        
# Return the mask

    return mask
        

#==============================================================================
# Test case to check if things are working
#==============================================================================
        
def test():

    from nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist
    
    dir_out    = '/Users/throop/Data/NH_Jring/out/'

    objectfile = 'lor_0034604523_0x630_sci_1_opnav_objects.txt'
    objectfile_base = objectfile.split('/')[-1]
    path_objectfile = dir_out + objectfile_base
    
    file_image = objectfile_base.replace('_objects', '').replace('.txt', '.fit')
    dir_images = '/Users/throop/Data/NH_Jring/data/jupiter/level2/lor/all/'
    dir_out    = '/Users/throop/Data/NH_Jring/out/'

    im = hbt.read_lorri(dir_images + file_image)

# Create the mask

    mask = nh_jring_mask_from_objectlist(objectfile)    

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile.     
    
    im_masked = im.copy()
    im_masked[mask == True] = np.median(im)
    
    hbt.figsize((10,10))

# And make some plot to show if it matches up, or not

    plt.imshow(stretch(im_masked))
    plt.show()
    plt.imshow(stretch(im)) 
    plt.imshow(stretch(hbt.lorri_destripe(im)))
    plt.show()
