#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:21:34 2016

@author: throop
"""

# Now that we have generated the radial profiles in nh_jring_gui.py, read and plot the results.

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
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astroquery.vo_conesearch import conesearch
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import time
from scipy.interpolate import griddata

import re # Regexp
import pickle # For load/save


# HBT imports

import hbt

hbt.figsize((12,8))
plt.set_cmap('Greys_r')

dir_input = '/Users/throop/data/NH_Jring/out/'
files_input = glob.glob(dir_input + '*_export.pkl')

num_radius         = 300
num_azimuth        = 300
    
numfiles           = np.size(files_input)
ang_phase          = np.zeros((numfiles))
ang_elev           = np.zeros((numfiles))
ew                 = np.zeros((numfiles))
profile_radius_dn  = np.zeros((numfiles, num_radius))
profile_azimuth_dn = np.zeros((numfiles, num_azimuth))
et                 = np.zeros((numfiles))
index_image        = np.zeros((numfiles)).astype(int)
index_group        = np.zeros((numfiles)).astype(int)
azimuth            = np.zeros((numfiles, num_azimuth))
radius             = np.zeros((numfiles, num_radius))

image              = np.zeros((numfiles, num_radius, num_azimuth))

dist_inner = 128000
dist_outer = 132000

#==============================================================================
# Read in all the files.
# Extract and process the images, radial profiles, phase angle, et, etc.
#==============================================================================
for i,file in enumerate(files_input):
    lun = open(file, 'rb')
    vals = pickle.load(lun)
    lun.close

    (image[i,:,:], et[i], radius[i,:], azimuth[i,:], profile_radius_dn[i,:], profile_azimuth_dn[i,:], \
       ang_elev[i], ang_phase[i], 
       index_image[i], index_group[i]) = vals

#==============================================================================
# Do some one-off radial extractions, for cases that need some customization
#==============================================================================

for i in range(numfiles):

    if (index_group[i] in [5,6,7]):
        hbt.figsize((10,5))
        
        # Set limits for where the extraction should happen
        
        # For the radial profile, we take only a portion of the unwrapped frame. 
        # e.g., the inner 30%. We exclude the region on the edges, since it is smeared, and has contamination.
        # We are not starved for photons. Take the best portion of the signal and use it.
        
        frac_profile_radial = np.array([0.35, 0.65])
        
        if (index_group[i] == 5):
            frac_profile_radial = np.array([0.1, 0.4]) # What fractional region of azimuthal range do we use for radial profile
#            frac_profile_radial = np.array([0.7, 0.95]) # What fractional region of azimuthal range do we use for radial profile
            
        if (index_group[i] == 7) and (index_image[i] == 37):
            frac_profile_radial = np.array([0.05, 0.25])
                           
        # For azimuthal profile, focus on the main ring. Don't focus on the diffuse inner region.
        # It is harder to do that photometry, more artifacts, fainter, and probalby more spread out anyhow.
        
        # Define distances of [outer box, outer ring, inner ring, inner box] in km
        # These four radii define the position of the inner and outer regions to subtract as bg, 
        # and the inner region to use for the signal.
        
        limits_profile_azimuth = np.array([131e3,129.5e3,127e3,126e3]) # Radial limits, used for az profile
        
        limits_profile_azimuth_bins = limits_profile_azimuth.astype(int).copy() * 0

        # XXX I think ths should be rewritten with hbt.x2bin()
        
        for j,r in enumerate(limits_profile_azimuth):
            limits_profile_azimuth_bins[j] = int(hbt.wheremin(abs(radius[i,:] - r)))
        
        limits_profile_radial_bins = np.shape(image[i,:,:])[1] * frac_profile_radial

        limits_profile_radial_bins = limits_profile_radial_bins.astype(int)            

        #==============================================================================
        # Extract radial and azimuthal profiles, using subsections
        #==============================================================================
        
        # I am using the *mean* here along each row and column. 
        
        #profile_azimuth_bg_inner = np.nansum(dn_grid[limits_profile_azimuth_bins[1]:limits_profile_azimuth_bins[0],:],0)
                
        profile_azimuth_bg_inner = np.nanmean(image[i,limits_profile_azimuth_bins[1]:limits_profile_azimuth_bins[0],:],axis=0)
        profile_azimuth_core     = np.nanmean(image[i,limits_profile_azimuth_bins[2]:limits_profile_azimuth_bins[1],:],axis=0)
        profile_azimuth_bg_outer = np.nanmean(image[i,limits_profile_azimuth_bins[3]:limits_profile_azimuth_bins[2],:],axis=0)

        # Get profile in DN
        profile_azimuth_central  = profile_azimuth_core - (profile_azimuth_bg_inner + profile_azimuth_bg_outer)/2
        profile_radius_central   = np.nanmean(image[i,:,limits_profile_radial_bins[0]:limits_profile_radial_bins[1]],1)
        
        extent = [azimuth[i,0], azimuth[i,-1], radius[i,0], radius[i,-1]]

        f = (np.max(radius[i,:]) - np.min(radius[i,:])) / (np.max(azimuth[i,:]) - np.min(azimuth[i,:]))
        aspect = 0.5/f # Set the aspect ratio

#        stretch = astropy.visualization.PercentileInterval(self.stretch_percent)  # PI(90) scales array to 5th .. 95th %ile. 

        # For now, we are scaling this vertically by hand. Might have to revisit that.

        # Make a 2 x 1 plot: image on LHS, profile on RHS
        
        plt.subplot(1,2,1)
        plt.imshow(image[i,:,:], extent=extent, aspect=aspect, vmin=-15, vmax=20, origin='lower') # aspect='auto'a
        plt.title("{}/{}: {}".format(index_group[i], index_image[i], files_input[i].split('/')[-1]))
        plt.hlines(limits_profile_azimuth[0], -10, 10, color='purple')
        plt.hlines(limits_profile_azimuth[1], -10, 10, color='purple')
        plt.hlines(limits_profile_azimuth[2], -10, 10, color='purple')
        plt.hlines(limits_profile_azimuth[3], -10, 10, color='purple')
#        self.ax3.set_xlim(hbt.mm(bins_azimuth))
        
        plt.vlines(azimuth[i,limits_profile_radial_bins[0]],-1e10, 1e10)
        plt.vlines(azimuth[i,limits_profile_radial_bins[1]],-1e10, 1e10)

        plt.ylim((radius[i,0], radius[i,-1]))
        plt.xlim((azimuth[i,0], azimuth[i,-1]))
        
        plt.xlabel('Azimuth [radians]')
        plt.ylabel('Radius [$R_J$]')
        
        plt.subplot(1,2,2) # Plot the radial profile on RHS
        plt.plot(radius[i,:], profile_radius_central)
        plt.xlim((122000,133000))
        plt.show()
        
        # And now export it to the main array
        
        profile_radius_dn[i,:] = profile_radius_central
        profile_azimuth_dn[i,:] = profile_azimuth_central
        
#==============================================================================
# Cross-correlate these signals and shift them radially in order to align them using radial profiles
#==============================================================================

limit_shift = np.array([127000,130000])
limit_shift_bin = hbt.x2bin(limit_shift, radius[0,:])

shift = np.zeros((numfiles)).astype(int)
profile_radius_dn_roll = profile_radius_dn.copy()

dx_max = 20  # look at all the posible shifts, from -dx_max to +dx_max, which is 2*dx_max+1 different integers

correl = np.zeros((numfiles, dx_max*2 + 1)) 

for i in range(numfiles):
    for j,dx in enumerate(hbt.frange(-dx_max, dx_max).astype(int)):
        correl[i,j] = np.correlate(profile_radius_dn[0,limit_shift_bin[0]:limit_shift_bin[1]], \
                           np.roll(profile_radius_dn[i,limit_shift_bin[0]:limit_shift_bin[1]],dx))
        if np.isnan(correl[i,j]): # Some of the radial profiles have nans in them. For these, just skip the correlation 
            correl[i,j] = 1

    shift[i] = (hbt.wheremax(correl[i,:])) - dx_max
    profile_radius_dn_roll[i,:] = np.roll(profile_radius_dn[i,:],shift[i])

# Shift each image vertically

for i in range(numfiles):
    bin_inner_vnorm = hbt.x2bin(131000, radius[i,:])
    bin_outer_vnorm = hbt.x2bin(133000, radius[i,:])
    profile_radius_dn_roll[i,:] -= np.mean(profile_radius_dn_roll[i,bin_inner_vnorm:bin_outer_vnorm])


# Fit and subtract a polynomial from each curve
    
# Make a plot of all the radial profiles, now aligned both vertically and horizontally

dy = 3

hbt.figsize((4,8))

for i in range(numfiles):    
    plt.plot(radius[i,:], profile_radius_dn_roll[i,:] + i * dy)
    
plt.xlim((125000, 130000))
plt.show()


for i in range(numfiles):    
    plt.plot(radius[i,:], profile_radius_dn[i,:] + i * dy)
    
plt.xlim((125000, 130000))
plt.show()


xlim = (126000,131000)
xlim_bin = hbt.x2bin(xlim, radius[0,:])

xlim_ew = (127700,129200)
xlim_ew_bin = hbt.x2bin(xlim_ew, radius[0,xlim_bin[0]:xlim_bin[1]])


for i in range(numfiles):    
    arr = profile_radius_dn_roll[i,xlim_bin[0]:xlim_bin[1]]
    try:
        arr_fit = hbt.remove_polyfit(arr, degree=5)
    except (ValueError):
        arr_fit = 0 * arr
        
    plt.plot(radius[i,xlim_bin[0]:xlim_bin[1]], arr_fit + i * dy)

    ew[i] = np.sum(arr_fit[xlim_ew_bin[0] : xlim_ew_bin[1]])
    
plt.xlim((125000, 130000))
plt.show()


## Calculate EW
#
#for i in range(numfiles):
#    bin_inner = hbt.x2bin(129000, radius[i,:])
#    bin_outer = hbt.x2bin(131000, radius[i,:])
#    ew[i] = np.sum(profile_radius_dn_roll[i,bin_inner:bin_outer])
    
# Make a plot of EW vs. phase

hbt.figsize((5,5))
plt.plot(ang_phase*hbt.r2d, ew, marker = 'o', linestyle='none')
plt.xlabel('Phase Angle [deg]') 
plt.ylabel('EW [DN * km]')
plt.show()

# Same plot, but linear-log

hbt.figsize((7,5))
plt.plot(ang_phase*hbt.r2d, ew, marker = 'o', linestyle='none')
plt.xlabel('Phase Angle [deg]') 
plt.ylabel('EW [DN * km]')
plt.yscale('log')
plt.xlim((0,180))
plt.ylim((1e-1,1e3))
plt.show()

# Make a plot of all of the images, individually, along with radial profiles

plt.set_cmap('Greys_r')
stretch = astropy.visualization.PercentileInterval(1)  # PI(90) scales array to 5th .. 95th %ile. 

#hbt.figsize((1,20)) # Horizontal, vertical

hbt.figsize((10,3*numfiles))      # x, y
plt.show()

for i in range(numfiles):
    
    extent = [azimuth[i,0], azimuth[i,-1], radius[i,0], radius[i,-1]]
    f = (np.max(radius[i,:]) - np.min(radius[i,:])) / (np.max(azimuth[i,:]) - np.min(azimuth[i,:]))
    
    aspect = 0.4/f # Set the aspect ratio
    
    plt.subplot(numfiles,2,2*i+1)  # y, x, n
    plt.imshow(image[i,:,:], aspect=aspect, vmin=-10,vmax=20,extent=extent, origin='lower')
    plt.title("i={}: {}/{}: {}".format(i, index_group[i], index_image[i], \
                                       files_input[i].split('/')[-1].replace('_export.pkl', '')))

    plt.subplot(numfiles,2,2*i+2)
    plt.plot(radius[i,:], profile_radius_dn_roll[i,:])
    plt.xlim((122000,130000))
    
plt.show()    

#==============================================================================
# Plot a few individual profiles where we can see this three-humped ring structure
#==============================================================================

hbt.figsize((5,3))
groupnum = np.array([5,8,8,8])
imagenum = np.array([2,55,70,106])

j=0

for i in ((0,1,4,15,16,33)):
#    mask = (groupnum[i] == index_group) & (imagenum[i] == index_image)   # '&' works, but 'and' does not
#    if (np.sum(mask) > 0):
        index = np.where(mask == True)[0][0]
        plt.plot( radius[i,:], profile_radius_dn_roll[i,:] + j)
        j+=1
        
#        plt.imshow(image[index], origin='lower')
#        plt.ylim((150,250))
#        plt.title("{}/{}".format(index_group[index], index_image[index]))

plt.xlim((125000,130000))
plt.show()
r2d = hbt.r2d

fs = 12

hbt.figsize((10,6))
plt.plot( radius[0,:], profile_radius_dn_roll[0,:] + j, label = repr(ang_phase[0]*r2d))
#plt.plot( radius[0,:], profile_radius_dn_roll[3,:] + j)
#plt.plot( radius[0,:], profile_radius_dn_roll[9,:] + j)
plt.plot( radius[0,:], profile_radius_dn_roll[15,:] + j, label = repr(ang_phase[15]*r2d))
plt.plot( radius[0,:], profile_radius_dn_roll[20,:] + j, label = repr(ang_phase[20]*r2d))
plt.plot( radius[0,:], profile_radius_dn_roll[28,:] + j, label = repr(ang_phase[28]*r2d))
plt.plot( radius[0,:], profile_radius_dn_roll[33,:] + j, label = repr(ang_phase[33]*r2d))
plt.xlim((125000,130000))
plt.ylabel('DN', fontsize=fs)
plt.ylim((0,36))
#plt.legend() # 5, 20 93, 140
plt.xlabel('Orbital Distance [km]')

plt.show()

        