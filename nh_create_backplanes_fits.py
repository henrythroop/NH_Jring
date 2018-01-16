#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

File to create the backplanes for MU69 encounter images.
Written for hazard search ORTs, but should work more generally.

Created on Tue Jan  9 13:19:09 2018

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

from astropy.coordinates import SkyCoord

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt
from   image_stack import image_stack
from   create_backplane import create_backplane

def nh_create_backplanes_fits(file_in = None, 
                              file_out = None, 
                              name_target = 'MU69', 
                              name_observer = 'New Horizons', 
                              clobber = False,
                              do_plot = True,
                              type = 'Sunflower'):
    
    """
    Function to create backplanes and add them to a FITS file.
    
    Idea is that this should be runnable as a commandline script.
    
    Parameters
    ----
    
    file_in:
        Input FITS file. Should be fully navigated and have good header info (correct WCS, etc)
    file_out:
        Output FITS file to be written. If None, then no file is written.
    name_target:
        String. Name of central body of the backplanes. For now, must be 'MU69'.
    name_observer:
        String. Name of observer. Usually 'New Horizons'
    clobber:
        Boolean. The original file will be overwritten, iff file_out == file_in and file==True.
    type:
        Type of orbit to assume for the backplanes. Can be 'Sunflower'. Other orbits might be added later, like ones 
        that face north pole, or are tilted relative to ecliptic, rather than face the Sun.
    
    Output
    ----    
    
    tuple:
        A tuple of(radius, azimuth, d_ra, d_dec, etc). This has all of the backplanes in it.
    
    Return status
    ----
        0 or 1. 0 = no problems.
        
    """

#- Open the file
#- Read the time, mode (1x1 vs 4x4), image scale, etc.
#- Compute the backplanes
#- Copy the FITS fileconso
#- Write the new FITS file
#- Return output to user

# Intialize graphics
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    plt.set_cmap('Greys_r')

    frame = 'J2000'
    name_target = 'MU69'
    name_observer = 'New Horizons'        
    abcorr = 'LT' # Not sure which abcorr to use, but LT vs LT+S makes ~ 1-pixel difference
    
# Start up SPICE
    
    file_kernel = 'kernels_kem.tm'
    sp.furnsh(file_kernel)
    
# Load the image
        
    hdu = fits.open(file_in) 
    data = hdu['PRIMARY'].data
    header = hdu['PRIMARY'].header

# Grab some critical info from header

    et      = header['SPCSCET']
    sformat = header['SFORMAT']
    dx_pix  = header['NAXIS1']
    dy_pix  = header['NAXIS2']
    mode    = header['SFORMAT']
    exptime = header['EXPTIME'] 

# Do some processing on the header

    utc = sp.et2utc(et, 'C', 0)
    w = WCS(header)
    
# Look up position of MU69 in pixels.
    
    (st,lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
    vec_obs_mu69 = st[0:3]
    (_, ra, dec) = sp.recrad(vec_obs_mu69)
    (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)

    dist_mu69 = sp.vnorm(st[0:3])*u.km
    
# Display the image, and MU69.

    if do_plot:
        
        radec_str = SkyCoord(ra=ra*u.rad, dec=dec*u.rad).to_string('dms')
    
        hbt.figsize((10,10))
        hbt.set_fontsize(12)
        plt.subplot(projection=w)    
        plt.imshow(stretch(data))
        plt.plot(pos_pix_x, pos_pix_y, color='orange', ls='None', marker = 'o', ms=7, alpha=1, 
                 label = 
           '{}, RA={:0.3f}, Dec={:0.3f}, {}'.format(name_target, ra*hbt.r2d, dec*hbt.r2d, radec_str))
           
    # Plot a sunflower ring, using the Sunflower rotating reference frame
    
        frame_sunflower = '2014_MU69_SUNFLOWER_ROT'  # INERT frame is not working yet. 
    
        az_ring = hbt.frange(0,2*math.pi, 100) # Set up azimuth angle
        radius_ring = 10000*u.km
     
        (n,re) = sp.bodvrd( 'MU69', 'RADII', 3)
        mx_sunflower_j2k = sp.pxform(frame_sunflower, 'J2000', et)
        radius = radius_ring.to('km').value
    
    # For each point in the ring, just define an XYZ point in space.
    # Then, we will call PXFORM to transform from Sunflower frame, to J2K frame.
    # As per MRS sketches, sunflower frame is defined s.t. ring is in the X-Z plane.
                                
        for az_ring_i in az_ring:
    
            vec_ring_i_mu69 = [radius * np.cos(az_ring_i), 0, radius * np.sin(az_ring_i)] # Circle in X-Z plane. 
            vec_mu69_ring_i_j2k  = sp.mxvg(mx_sunflower_j2k, vec_ring_i_mu69, 3,3)
            vec_obs_ring_i_j2k = vec_mu69_ring_i_j2k + vec_obs_mu69
            (_, ra, dec) = sp.recrad(vec_obs_ring_i_j2k)
            (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
            plt.plot(pos_pix_x, pos_pix_y, color='orange', ls='None', marker = 'o', ms=1, alpha=1)
        
        plt.title('{}, {}, {}'.format(os.path.basename(file_in), utc, abcorr))
        plt.legend(framealpha=1)
        plt.show()
        
        print("MU69 location from SPICE: RA {:0.3f}, Dec {:0.3f} → {}".format(ra*hbt.r2d, dec*hbt.r2d, radec_str))

# =============================================================================
# Call a routine to actually create the backplanes, which returns them as a tuple.
# =============================================================================

    (planes, descs) = create_backplane(file_in,
                              type = 'Sunflower',
                              name_target = name_target,
                              name_observer = name_observer
                              )
  
# =============================================================================
# Now write everything to a new FITS file. 
# =============================================================================
    
    file_out = file_in.replace('.fit', '_backplaned.fit')   # Works for both .fit and .fits
    
    # Open the existing file
    
    hdu = fits.open(file_in)
    
    # Go thru all of the new backplanes, and add them one by one. For each, create an ImageHDU, and then add it.
    
    for key in planes.keys():
        hdu_new = fits.ImageHDU(data=planes[key].astype(np.float32), name=key, header=None)
        hdu.append(hdu_new)
    
    # Add relevant header info
    
    hdu[0].header['COMMENT'] = '*********************************************************'
    hdu[0].header['COMMENT'] = '*** BACKPLANE INFO                                    ***'
    hdu[0].header['COMMENT'] = '*********************************************************'
    
    for i,desc in enumerate(descs):
        hdu[0].header['BKPLN_{}'.format(i)] = desc
        
    # Write to a new file
    
    hdu.writeto(file_out, overwrite=True)

    print("Wrote: {}; {} planes; {:.1f} MB".format(file_out, 
                                                   len(hdu), 
                                                   os.path.getsize(file_out)/1e6))

    hdu.close()
    
# Return to user
    
    return(0)

# =============================================================================
# Do some tests, to validate the function
# =============================================================================

if (__name__ == '__main__'):
    
    file_in = os.path.join(os.path.expanduser('~'), 'Data', 'NH_KEM_Hazard', 'ORT1_Jan18', 
#                                'lor_0406731132_0x633_sci_HAZARD_test1.fit'
                                'lor_0406731132_0x633_sci_HAZARD_test1-1.fit'
#                                'lor_0368310087_0x633_sci_HAZARD_test1.fit'
                                )

    
    # Create the backplanes
    
    nh_create_backplanes_fits(file_in, None, do_plot=False)
