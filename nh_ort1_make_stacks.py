#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 22:36:55 2018

@author: throop
"""

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
import warnings
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
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)
import imreg_dft as ird                    # Image translation

from astropy.coordinates import SkyCoord

import os.path

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   get_radial_profile_circular import get_radial_profile_circular 

def nh_ort1_make_stacks():
    
    """
    This program takes a directory full of individual NH KEM Hazard frames, stacks them, and subtracts
    a stack of background field. This reveals rings, etc. in the area.
    
    Written for NH MU69 ORT1, Jan-2018.  
    """
    
    do_force = False  # Boolean: Do we force reloading of all of the images from FITS, or just restore from pkl?
                      # Pkl (aka False) is faster. But if we have made changes to the core algorithms, must
                      # reload from disk (aka True).
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
    reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
    
    dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'

    zoom = 4
    
    # Set the edge padding large enough s.t. all output stacks will be the same size.
    # This value is easy to compute: loop over all stacks, and take max of stack.calc_padding()[0]
    
    padding = 61   
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
        
    # Set the RA/Dec of MU69. We could look this up from SPICE but it changes slowly, so just keep it fixed for now.
    
    radec_mu69 = (4.794979838984583, -0.3641418801015417)
    
    # Load and stack the field images
    
    stack_field = image_stack(os.path.join(dir_data, reqid_field), do_force=do_force)
    stack_field.align(method = 'wcs', center = radec_mu69)
    img_field  = stack_field.flatten(zoom=zoom, padding=padding)
    
    if do_force:
        stack_field.save()

    hbt.figsize((12,12))
    hbt.set_fontsize(15)
    
    for reqid in reqids_haz:
        stack_haz = image_stack(os.path.join(dir_data, reqid), do_force=do_force)
        stack_haz.align(method = 'wcs', center = radec_mu69)
        img_haz  = stack_haz.flatten(zoom=zoom, padding=padding)

        if do_force:
            stack_haz.save()
            
        # Make the plot
        
        diff = img_haz - img_field
        diff_trim = hbt.trim_image(diff)
        plt.imshow(stretch(diff_trim))
        plt.title(f"{reqid} - field, zoom = {zoom}")

        # Save the stacked image as a FITS file
        
        file_out = os.path.join(dir_data, reqid, "stack_n{}_z{}.fits".format(stack_haz.size[0], zoom))
        hdu = fits.PrimaryHDU(stretch(diff_trim))
        hdu.writeto(file_out, overwrite=True)
        print(f'Wrote: {file_out}')    
        
        # Save the stack as a PNG
        
        file_out_plot_stack = file_out.replace('.fits', '.png')
        plt.savefig(file_out_plot_stack, bbox_inches='tight')
        print("Wrote: {}".format(file_out_plot_stack))

        # Display it 
        # This must be done *after* the plt.savefig()
        
        plt.show()
        
        # Make a radial profile
        
        pos =  np.array(np.shape(diff))/2
        (radius, profile) = get_radial_profile_circular(diff, pos=pos, width=1)
    
        hbt.figsize((10,8))
        hbt.set_fontsize(15)
        plt.plot(radius, profile)
        plt.xlim((0, 50*zoom))
        plt.ylim((-1,np.amax(profile)))
        plt.xlabel('Radius [pixels]')
        plt.title(f'Ring Radial Profile, {reqid}, zoom={zoom}')
        plt.ylabel('Median DN')
        plt.show()

# =============================================================================
# Calculate how many DN MU69 should be at encounter (K-20d, etc.)
# Or alternatively, convert all of my DN values, to I/F values
# =============================================================================

        # Convert DN values in array, to I/F values
            
        RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
        RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
        
        C = profile # Get the DN values of the ring. Typical value is 1 DN.
        
        # Define the solar flux, from Hal's paper.
        
        FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
        F_solar = FSOLAR_LORRI # Flux from Hal's paper
        
        RSOLAR = RSOLAR_LORRI_4X4
        
        # Calculate the MU69-Sun distance, in AU (or look it up). 
        
        km2au = 1 / (u.au/u.km).to('1')
        
        et = stack_haz.t['et'][0]
        (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
        r_nh_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
        
        (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'Sun')
        r_sun_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
        
        pixscale_km =  (r_nh_mu69/km2au * (0.3*hbt.d2r / 256)) / zoom # km per pix (assuming 4x4)
        
        TEXP = stack_haz.t['exptime'][0]
        
        I = C / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.
        
        # Apply Hal's conversion formula from p. 7, to compute I/F and print it.
        
        IoF = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper
        
        plt.plot(radius * pixscale_km, IoF)
        plt.xlim((0, 50000))
        plt.ylim((-1e-7, 4e-7))
#        plt.ylim((0,np.amax(IoF)))
#        plt.yscale('log')
        plt.xlabel('Radius [km]')
        plt.title(f'Ring Radial Profile, {reqid}, zoom={zoom}')
        plt.ylabel('Median I/F')
        file_out_plot_profile = file_out.replace('.fits', '_profile.png')
        plt.savefig(file_out_plot_profile, bbox_inches='tight')
        plt.show()
        print(f'Wrote: {file_out_plot_profile}')
        
        # Write it to a table
        t = Table([radius, radius * pixscale_km, profile, IoF], names = ['RadiusPixels', 'RadiusKM', 'DN/pix', 'I/F'])
        file_out_table = file_out.replace('.fits', '_profile.txt')
        t.write(file_out_table, format='ascii', overwrite=True)
        print("Wrote: {}".format(file_out_table))
        
#    # Now convert to 'normal I/F'
#    
#    # Define mu = cos(e), where e = emission angle, and e=0 is face-on.
#    
#    e = 0 # Assume ring is face-on. The sunflower orbits are. 
#    
#    mu = np.cos(e)
#    
#    #        mu_2D = np.transpose(np.tile(mu, (self.num_bins_radius(),1)))
#    
#    # Calculate the normal I/F
#    
#    IoF_normal = 4 * mu * IoF
        
if (__name__ == '__main__'):
    nh_ort1_make_stacks()
    
    
#a  = plt.imshow(hbt.dist_center(100))
#plt.savefig('test.png', bbox_inches='tight')
#plt.show()
    
    