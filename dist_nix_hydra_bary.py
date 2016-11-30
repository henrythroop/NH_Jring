#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 12:35:39 2016

This is a one-off program to make a simple plot of Nix / Hydra orbital distances vs. time.
This is used in my UV Ring Occultation code, where I want to know whether we crossed the orbits, or not.
Concl: it's complicated!

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

import hbt

file_tm = '/Users/throop/gv/dev/kernels/15sci_rhr_25aug16.tm' # Def need to use an encounter kernel here.

cspice.furnsh(file_tm)

day = 86400.
month = 30*day

num_dt = 1000
frame = 'J2000'
abcorr = 'LT'

utc_start = '2015 Jul 16 21:25:00'  # Start of OC3 stellar occ
et_start = cspice.utc2et(utc_start)
et_end   = et_start + 1*month

et = hbt.frange(et_start, et_end, num_dt)

d_hyd_bar = []  # distance
d_nix_bar = [] # distance nix - pluto barycenter
d_hyd_plu = [] # distance
d_nix_plu = []

for et_i in et:
    (st,lt) = cspice.spkezr('Nix', et_i, frame, abcorr, 'Pluto')
    d_nix_plu.append(cspice.vnorm(st[0:3]))
 
    (st,lt) = cspice.spkezr('Hydra', et_i, frame, abcorr, 'Pluto')
    d_hyd_plu.append(cspice.vnorm(st[0:3]))
    
    (st,lt) = cspice.spkezr('Nix', et_i, frame, abcorr, 'Pluto Barycenter')
    d_nix_bar.append(cspice.vnorm(st[0:3]))
    
    (st,lt) = cspice.spkezr('Hydra', et_i, frame, abcorr, 'Pluto Barycenter')
    d_hyd_bar.append(cspice.vnorm(st[0:3]))

    
offset = 10000
    
d_hyd_plu = np.array(d_hyd_plu)
d_nix_plu = np.array(d_nix_plu)
d_hyd_bar = np.array(d_hyd_bar)
d_nix_bar = np.array(d_nix_bar)

t_day = (et-et_start)/day
plt.plot(t_day, d_hyd_plu, label = 'Hydra from Pluto')
plt.plot(t_day, d_hyd_bar)
plt.plot(t_day, d_nix_plu + offset)
plt.plot(t_day, d_nix_bar + offset, label = 'Nix from Barycenter')
plt.xlabel('Days')
plt.ylabel('Distance [km]')
plt.legend()
plt.show()
    
    