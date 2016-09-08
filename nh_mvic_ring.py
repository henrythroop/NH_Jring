# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 22:55:34 2016

@author: throop
"""

dir = '/Users/throop/Data/NH_MVIC_Ring/'
file = 'mvic_d305_sum_mos_v1.fits'

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
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import astropy.visualization
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools   import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy     import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils   import daofind
import wcsaxes
import time
from scipy.interpolate import griddata

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

###

import hbt

dir  = '/Users/throop/Data/NH_MVIC_Ring' 

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

cspice.furnsh(file_tm) # Start up SPICE

#==============================================================================
# Look at the first file (small one)
#==============================================================================

files = ['mvic_d305_sum_mos_v1.fits', 'ringdep_mos_v1.fits']

file = dir + '/mvic_d305_sum_mos_v1.fits'
file = dir + '/mvic_d305_sum_mos_v1-new-image_fixed.fits' # Load the navigated image, with WCS

plt.rcParams['figure.figsize'] = 6,56
hdulist = fits.open(file)
im = hdulist['PRIMARY'].data

plt.set_cmap('Greys_r')
stretch = astropy.visualization.PercentileInterval(99)
plt.imshow(stretch(im))
#plt.imshow(im)

header = hdulist['PRIMARY'].header

#stop

#==============================================================================
# Look at the second file (larger)
#==============================================================================

#file = dir + files[1]
#
#plt.rcParams['figure.figsize'] = 16,156
#hdulist = fits.open(file)
#im = hdulist['PRIMARY'].data
#
#plt.set_cmap('Greys_r')
#stretch = astropy.visualization.PercentileInterval(99.8)
#plt.imshow(stretch(im))
#plt.show()

#==============================================================================
# Now generate a backplane
#==============================================================================

planes = hbt.create_backplane(file, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')

file_out = dir_out + file_short    
file_out = file_out.replace('.fits', '_planes.pkl')
    
#header['SPCSCET'] =     (490257526.5058528, "[s past J2000] Spacecraft mid-obs time, TDB")   

#file_out = file.replace('.fits', '_fixed.fits')
#hdulist.writeto(file_out)

    # Write one variable to a file        

file_out = file_out.replace('.fit', '_planes.pkl')

lun = open(file_out, 'wb')
pickle.dump(plane, lun)
lun.close()

print "Wrote: " + file_out
print
    