#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 13:56:13 2018

read_mehoke_sizes

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

from   astropy.utils import data
from   astropy.wcs import WCS

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
from   astroquery.vo_conesearch import conesearch
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

from   matplotlib.figure import Figure
from   matplotlib.patches import Rectangle

from   pymiecoated import Mie

# HBT imports

import hbt

def read_mehoke_sizes():
    
    """
    This file reads a .csv file which I created based on Doug Mehoke's Excel file
    "Damge Areas.xls." [sic] His file lists the s/c damage susceptibility for New Horizons.
    I process this data, sort it, and make some plots. This is used in my back-of-envelope
    hazard plots -- e.g., created by nh_ringhazard.py .
    """
    
    hbt.figsize((10,8))
    
    file_in = 'DamageAreasMehoke.csv'  # sic: Damge not Damage
    
    t = Table.read(file_in, format = 'csv')
    t.sort('r_c um')
    
    plt.plot(t['Area in2'], t['r_c um'])
    
    area = np.array(t['Area in2']) * (2.54)**2 * (u.cm)**2  # Convert to cm2
    r_c = np.array(t['r_c um']) * u.micron
    component = t['Component']
    
    fac_projection = 0.35    # Roughly, what fraction of the s/c is pointed in the ram 
                             # direction at any time. This is very rough. It is somewhere
                             # between 1/6 and 1/2.
                             
    area_cum = np.cumsum(area)
    
    plt.plot(r_c, area, linestyle = 'none', marker = 'o')
    plt.xlabel('r [micron]')
    plt.ylabel('Component Area, cm2')
    plt.show()
    
    plt.plot(r_c, area_cum.to(u.m**2) * fac_projection, linestyle = '-', marker = 'o', ms=3)
    plt.xlabel(r'Particle Radius [$\mu m$]')
    plt.ylabel(r'Projected S/C Area Susceptible to Damage by < r [$m^2$]')
    plt.ylim((1e-3, 5))
    plt.yscale('log')
    plt.show()

if (__name__ == '__main__'):
    read_mehoke_sizes()
    