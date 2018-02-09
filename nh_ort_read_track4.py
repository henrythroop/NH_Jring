#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 16:24:19 2018

Program to read the ORT 'Track 4' results -- that is, the ring flythru simulations of 
HBT + MRS.

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

from   datetime import datetime

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

dir_mrs = '/Users/throop/Data/ORT1/showalter'

files_mrs = glob.glob(os.path.join(dir_mrs, 'ort1*dust'))

for file in files_mrs:
    print(file)
    t = Table.read(file, format='ascii', names = ('delta_et', 'n_0', 'n_1', 'n_2', 'n_3', 'n_4', 'n_5', 'n_6'))

    for col in (t.colnames)[1:]:
        plt.plot(t['delta_et'], t[col], label = col)
        
    plt.yscale('log')
    plt.ylim(1,1e6)
    plt.xlim((-100,100))
    plt.legend(loc = 'upper left')
    plt.show()

    