#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:10:52 2018

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

# =============================================================================
# This is just a Q&D short routine to read in a directory of images, and plot them.
# I wrote this for the MU69 NYT ORT, but it coudl be used for anything.
# =============================================================================

# Start up SPICE if needed

hbt.figsize((10,10))
if (sp.ktotal('ALL') == 0):
    sp.furnsh('kernels_kem_prime.tm')
        
stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

dir_images = '/Users/throop/Data/ORT_Sep18/day2/lor/'

files = glob.glob(os.path.join(dir_images, '*'))

do_transpose = False

files = ['/Users/throop/Data/NH_Jring/data/jupiter/level2/lor/all/lor_0034715072_0x630_sci_1.fit']

for file in files:
    lun = fits.open(file)
    im = lun['PRIMARY'].data
    et = lun['PRIMARY'].header['SPCSCET']
    utc = sp.et2utc(et, 'C', 0)
    if do_transpose:
        im = np.transpose(im)
    plt.imshow(stretch(im), origin='lower')
    file_short = file.split('/')[-1]
    plt.title(f'{file_short} {utc}')
    plt.show()
    