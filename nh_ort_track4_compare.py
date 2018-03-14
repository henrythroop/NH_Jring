#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 16:24:19 2018

Program to read the ORT 'Track 4' results -- that is, the ring flythru simulations of 
HBT + MRS.

This routine is not used as part of the Hazard pipeline directly. Rather, it is used to compare
results of HBT and MRS, to make sure we are sending similar things to Doug Mehoke.

Track 4 are the lists of ring particle number densities as NH passes by MU69 on its path.
These files are created by HBT and MRS (independently). They are then given to Doug Mehoke,
for his use in Track 5.

@author: throop
"""

import pdb
import glob
import math
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

class track4_profile:
    
    def __init__(self, file):
        
            # Get the speed relative to MU69

        file_tm = 'kernels_kem_prime.tm'
        
        # Start up SPICE if needed
        
        if (sp.ktotal('ALL') == 0):
            sp.furnsh(file_tm)    
            
        utc_ca = '2019 1 Jan 05:33:00'
        et_ca  = sp.utc2et(utc_ca) 
        (st,lt) = sp.spkezr('New Horizons', et_ca, 'J2000', 'LT', 'MU69')
        
        velocity = sp.vnorm(st[3:6])*u.km/u.s
    
        # Save the velocity (relative to MU69)

        self.velocity = velocity

        # Save the name of file to read
        
        self.file = file

        # Save the area of the s/c
        
        self.area_sc = (1*u.m)**2
        
        return

# =============================================================================
# Read one of Mark Showalter's profiles. This is very easy.
# =============================================================================
    
    def read_mrs(self):
            
        t = Table.read(self.file, format='ascii', 
                       names = ('delta_et', 'n_0', 'n_1', 'n_2', 'n_3', 'n_4', 'n_5', 'n_6'))
    
        # Read the first row of the file as the size bins, and remove it.
            
        r_dust  = [t[0][1], t[0][2], t[0][3], t[0][4], t[0][5], t[0][6], t[0][7]] * u.mm
        t.remove_row(0)
        
        t['n_0'].unit = '1/km**3'
        t['n_1'].unit = '1/km**3'
        t['n_2'].unit = '1/km**3'
        t['n_3'].unit = '1/km**3'
        t['n_4'].unit = '1/km**3'
        t['n_5'].unit = '1/km**3'
        t['n_6'].unit = '1/km**3'
        
        n_dust = np.zeros((7))

        self.delta_et = t['delta_et']
        
        dt = (self.delta_et[1] - self.delta_et[0])*u.s   # Time incrememnt, in seconds

        t.remove_column('delta_et')

        for i in range(len(n_dust)):
            n_dust[i] = t[f'n_{i}'].sum()
        
        self.r_dust = r_dust
        self.t  = t
        self.dt = dt
        self.n_dust = n_dust
        
        self.process()
        print(f"Finished reading file MRS {self.file}")

# =============================================================================
# Read HBT profiles (new naming convention, for ORT2)
# =============================================================================

    def read_hbt(self):
        
        t = Table.read(self.file, format='ascii', 
                       names = ('delta_et', 'n_0', 'n_1', 'n_2', 'n_3', 'n_4', 'n_5', 'n_6'))
            
        t['n_0'].unit = '1/km**3'
        t['n_1'].unit = '1/km**3'
        t['n_2'].unit = '1/km**3'
        t['n_3'].unit = '1/km**3'
        t['n_4'].unit = '1/km**3'
        t['n_5'].unit = '1/km**3'
        t['n_6'].unit = '1/km**3'
        
        r_dust  = [t[0][1], t[0][2], t[0][3], t[0][4], t[0][5], t[0][6], t[0][7]] * u.mm
        t.remove_row(0)
        dt = (t['delta_et'][1] - t['delta_et'][0])*u.s  # Time incrememnt, in seconds

        self.t   = t
        self.delta_et = t['delta_et']
        t.remove_column('delta_et')
        
        self.dt = dt
        self.r_dust = r_dust
        self.n_dust = 0

        self.process()
        print(f"Finished reading file HBT {self.file}")
            
# =============================================================================
# Read HBT profiles (old naming convention, for ORT1).
# =============================================================================

#    def read_hbt_old(self, file):
#        
#        t = Table.read(file, format='ascii', 
#                       names = ('delta_et', 'n_0', 'n_1', 'n_2', 'n_3', 'n_4', 'n_5', 'n_6'))
#    
##        r_dust_row = t[0]
#        
#        t['n_0'].unit = '1/km**3'
#        t['n_1'].unit = '1/km**3'
#        t['n_2'].unit = '1/km**3'
#        t['n_3'].unit = '1/km**3'
#        t['n_4'].unit = '1/km**3'
#        t['n_5'].unit = '1/km**3'
#        t['n_6'].unit = '1/km**3'
#        
#        r_dust  = [t[0][1], t[0][2], t[0][3], t[0][4], t[0][5], t[0][6], t[0][7]] * u.mm
#        t.remove_row(0)
#        dt = (t['delta_et'][1] - t['delta_et'][0])*u.s  # Time incrememnt, in seconds
#
#        self.t   = t
#        self.delta_et = t['delta_et']
#        t.remove_column('delta_et')
#        
#        self.dt = dt
#        self.r_dust = r_dust
#        self.n_dust = 0
#
#        self.process()
#        print(f"Finished reading file HBT {self.file}")
#        
# =============================================================================
# Process some numbers. 
# =============================================================================
        
    def process(self):

        # Put the number density in a good place
        
        # Shape (7, 7201)
        number_density = \
                         np.array([self.t['n_0'].data, 
                                   self.t['n_1'].data, 
                                   self.t['n_2'].data, 
                                   self.t['n_3'].data, 
                                   self.t['n_4'].data, 
                                   self.t['n_5'].data, 
                                   self.t['n_6'].data])
        
#        self.t.copy()
        
        # Calculate cumulative

        number_cum = number_density.copy()
        
        number_cum *= ((1/u.km)**3 * self.area_sc * self.velocity * self.dt).to('1').value
        
        number_cum = np.cumsum(number_cum, axis=1)  # Units is number/sec
        
        self.number_cum     = number_cum
        self.number_density = number_density

        # Calculate a total line-of-sight density (ie, optical depth)
        
        tau = (np.sum( np.sum(self.number_cum, axis=1) * (math.pi * self.r_dust**2) ) / self.area_sc).to('').value
        
# =============================================================================
# Plot the profile
# =============================================================================
        
    def plot(self):
        
        t      = self.t
        n_dust = self.n_dust
        r_dust = self.r_dust
        
        hbt.figsize((13,9))
        

        # Plot number density
                  
        plt.subplot(2,2,1)    
        for i,col in enumerate( (t.colnames)[1:] ):
            plt.plot(self.delta_et, t[col], label = "{:.2f} mm".format(r_dust[i]))
        
        plt.yscale('log')
        plt.ylim(1e-6,1e8)
        plt.xlim((-300,300))
        plt.legend(loc = 'upper left')
        plt.title(os.path.basename(self.file))
        plt.xlabel('t from C/A [sec]')
        plt.ylabel('# particles/km3')
        
        # Plot cumulative number
        
        plt.subplot(2,2,2)    
        for i in range(len(self.r_dust)):
            plt.plot(self.delta_et, self.number_cum[i,:], 
                     label = "{:.2f}".format(r_dust[i]))
        
        n_dust = np.zeros((7))
        for i in range(len(n_dust)):
            n_dust[i] = t[f'n_{i}'].sum()
        plt.yscale('log')
        plt.ylim((1e-6,1e9))
        plt.xlim((-300,300))
        plt.legend(loc = 'upper left')
        plt.title(os.path.basename(self.file))
        plt.xlabel('t from C/A [sec]')
        plt.ylabel('Cum # particles hit on SC')
                   
        # Make a plot of the size distribution by itself
        
        plt.subplot(2,2,3)
        plt.plot(r_dust, n_dust, marker = 'o', linestyle = '--')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('r [mm]')
        plt.ylabel('N per km3')
        plt.ylim((1e1, 1e9))
        plt.xlim((1e-2, 1e1))
        plt.show()
    
# =============================================================================
# Run the method
# =============================================================================
 
if __name__ == '__main__':

    do_mrs = True
    do_hbt = True
    
    if do_mrs:    
        dir_mrs = '/Users/throop/Data/ORT1/showalter'
        files_mrs = glob.glob(os.path.join(dir_mrs, 'ort1*dust'))
        files_mrs = files_mrs[0:3]
        file_mrs = files_mrs[0]
        for file_mrs in files_mrs:
            track = track4_profile(file_mrs)
            track.read_mrs()
            self = track
            track.plot()

    if do_hbt:    
        dir_hbt = '/Users/throop/Data/ORT1/throop/track4/'
        dir_hbt = '/Users/throop/Data/ORT2/throop/track4/'
        files_hbt = glob.glob(os.path.join(dir_hbt, '*ort2*dust'))
        file_hbt = files_hbt[0]
        for file in files_hbt:
            track2 = track4_profile(file_hbt)
            track2.read_hbt()
            self = track2
            track2.plot()
    
    file_hbt = '/Users/throop/Data/ORT2/throop/track4/traj1_ort2-ring_speed1_q35_pv1_rho1_inc2.dust'
    file_mrs = '/Users/throop/Data/ORT1/showalter/ort1-ring_traj1_speed1_pv2_rho3_q35.dust'
