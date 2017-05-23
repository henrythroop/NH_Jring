#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:57:17 2017

This routine read the Hedman & Starck phase cuves in HS16 paper.
This is for plotting vs. the J-ring phase curve, in Lauer Pluto rings paper.
It is a one-off routine, and I just grab a screenshot to use in the paper.

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
import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling


#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

# HBT imports

import hbt

file_in = 'hs16_table6.txt'

#%%

t = Table.read(file_in,format='ascii',header_start=4, data_start=5)

t['Contrast'] *= 1e-6   # It is in units of 1e-6

hbt.figsize((12,8))
plt.plot(180-t['theta_ave'], t['Contrast'], marker='o', linestyle = 'none')
plt.yscale('log')
plt.ylim([1e-7,2e-3])
plt.xlim([0, 200])
plt.title('I/F, G ring, Hedman & Starck 2016')
plt.show()

#%%
