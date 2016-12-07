#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:28:06 2016

@author: throop
"""

# Read LORRI outbound Pluto rings mosaics and calculate I/F limit from them.
# This is not a rigorous reduction -- just want to get an answer so I can write an MT
# about using LORRI vs. MVIC for unresolved KEM observations. 
# HBT 7-Dec-2016

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
import spiceypy as sp # was cspice
import skimage
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
import sys  # For stdout.write, without newline
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

from astropy import units as u
from astropy.coordinates import SkyCoord

# HBT imports
import hbt