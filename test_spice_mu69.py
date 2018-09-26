#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 23:52:48 2018

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

# =============================================================================
# Just a q&d code to test SPICE and phase angle.
# I am comparing NH-MU69 Phase angles to those from STK.
# Results below agree with GV, but not STK (or rather, disagree w audit list). 
# I think it's an STK issue.
# For Kelsi Singer / Oliver White 22-Sep-2018.
# =============================================================================

# HBT imports

import hbt
file_tm = 'kernels_kem_prime.tm'

sp.furnsh(file_tm)


utc_limits_arr = ["2018 1 Dec 5:00", "2018 1 Dec 12:00"]
utc_limits_arr = ["2019 1 Jan 5:00", "2019 1 Jan 6:00"]

et_limits_arr = []

for utc in utc_limits_arr:
    et_limits_arr.append(sp.utc2et(utc))

num_et = 100


et_arr = hbt.frange(et_limits_arr[0], et_limits_arr[1], num_et)

phase_arr = []
utc_arr   = []

for et in et_arr:
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT+S', 'New Horizons')
    vec_sc_mu69 = st[0:3]

    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT+S', 'Sun')
    vec_sun_mu69 = st[0:3]
    
    ang_phase = sp.vsep(-vec_sc_mu69, -vec_sun_mu69)
    
    phase_arr.append(ang_phase)
    
    utc_arr.append(sp.et2utc(et, 'C', 0))
    
    print(f'Phase angle = {ang_phase*hbt.r2d} deg at {utc_arr[-1]}')

phase_arr = np.array(phase_arr)
et_arr = np.array(et_arr)

d_et = f'Days past {utc_limits_arr[0]}'
d_et = (et_arr - np.amin(et_arr))/86400

hbt.figsize(8,8)
hbt.fontsize(12)
plt.plot(d_et, np.array(phase_arr)*hbt.r2d)
plt.xlabel(f'Days past {utc_limits_arr[0]}')
plt.ylabel('Phase [deg]')
plt.title('NH-Sun Phase Angle')
plt.text(0, 11.64, utc_limits_arr[0])
plt.text(np.amax(d_et), 11.73, utc_limits_arr[1])

plt.show()

