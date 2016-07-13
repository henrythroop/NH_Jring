# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:58:26 2016

READ_ALICE_RING_OCC
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
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
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

# HBT imports

import hbt

####

#sequence = 'O_RING_Occ3'
sequence = 'O_Ring_Occ2'

dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali/all'

file_list = glob.glob(dir_images + '/*fit')

#file_list = file_list[0:10]
t_all = []
count_rate_all = []

for file in file_list:
    data = hbt.read_alice(file)
#    print "."
    startmet = data['header_spect']['STARTMET']
    dt = data['header_count_rate']['SAMPLINT']
    count_rate = data['count_rate']
    
    t = startmet + dt * hbt.frange(0, np.shape(count_rate)[0]-1)
    
    print "Read file " + file + ", MET = " + repr(startmet) + ', N = ' + repr(len(count_rate))
    
    t_all.append(t)
    count_rate_all.append(count_rate)

count_rate_flat = np.array([item for sublist in count_rate_all for item in sublist])
t_flat =          np.array([item for sublist in t_all for item in sublist])

# Remove the last few samples

count_rate_flat = count_rate_flat

#plt.plot(t_flat, count_rate_flat, marker = '.', linestyle='none')

binning = 1000
fs = 15
count_rate_s = hbt.smooth_boxcar(count_rate_flat, binning)[binning:-binning]
t_s = t_flat[binning:-binning]

plt.plot(t_s - t_s[0], count_rate_s, marker = '.', linestyle='none', ms=0.1)
plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
plt.ylabel('Count Rate', fontsize=fs)
plt.xlabel('t since start [sec]', fontsize=fs)


stop
   