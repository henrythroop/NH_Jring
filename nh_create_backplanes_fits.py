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

import os.path
import os


import astropy
from   astropy.io import fits
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np

import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs

from   astropy.wcs import WCS
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search

from   matplotlib.figure import Figure

# HBT imports

import hbt
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   plot_backplanes import plot_backplanes

def nh_create_backplanes_fits(file_in, name_target, frame, name_observer, file_out,
                              do_clobber = False,
                              do_verbose = False,
                              do_plot    = False):    
    """
    Function to create backplanes and add them to a FITS file.
    
    Idea is that this should be runnable as a commandline script.
    
    SPICE must already be running.
    
    Parameters
    ----
    
    file_in:
        Input FITS file. Should be fully navigated and have good header info (correct WCS, etc)
    name_target:
        String. Name of central body of the backplanes.
    frame:
        The SPICE frame to be used for the backplane. Usually IAU_JUPITER, IAU_PLUTO, etc.
    name_observer:
        String. Name of observer. Usually 'New Horizons'
    file_out:
        Output FITS file to be written. If None, then no file is written.
    do_clobber:
        Boolean. Overwrite the output file?
    do_verbose:
        Boolean. 
    
    Output
    ----    
    
    All output is written to specified file. Nothing output to user.
    
    Return status
    ----
        0 or 1. 0 = no problems.
        
    """


# Initialize the output file, and check if it exists
    
    if not(file_out):
        file_out = file_in.replace('.fit', '_backplaned.fit')   # Works for both .fit and .fits
        
    if os.path.exists(file_out) and not(do_clobber):
        raise(FileExistsError)
        
# Load the input image
        
    hdu    = fits.open(file_in) 
    header = hdu['PRIMARY'].header

# Grab info from header

    et      = header['SPCSCET']
    utc     = sp.et2utc(et, 'C', 0)
    w       = WCS(header)

# =============================================================================
# Call a routine to actually create the backplanes, which returns them as a tuple.
# =============================================================================

    (planes, descs) = compute_backplanes(file_in, name_target, frame, name_observer)
      
# =============================================================================
# Now write everything to a new FITS file. 
# =============================================================================
    
    # Open the existing file
    
    hdu = fits.open(file_in)
    
    if (do_verbose):
        print("Read: {}".format(file_in))
        
    # Go thru all of the new backplanes, and add them one by one. For each, create an ImageHDU, and then add it.
    
    for key in planes.keys():
        hdu_new = fits.ImageHDU(data=planes[key].astype(np.float32), name=key, header=None)
        hdu.append(hdu_new)
    
    # Add relevant header info
    
    keys = list(planes.keys())
    for i,desc in enumerate(descs):
        hdu[0].header['BKPLN_{}'.format(i)] = "{}: {}".format(keys[i], desc)
    
    hdu[0].header['BKPLNFRM'] = (frame,       'Name of SPICE frame used for backplanes')
    hdu[0].header['BKPLNTRG'] = (name_target, 'Name of SPICE target used for backplanes')

    # Add a comment / divider. Not working, not sure why.
    
#    hdu[0].header['COMMENT'] = '*********************************************************'
#    hdu[0].header['COMMENT'] = '*** BACKPLANE INFO                                    ***'
#    hdu[0].header['COMMENT'] = '*********************************************************'
#    
    # Write to a new file
    
    hdu.writeto(file_out, overwrite=True)

    if (do_verbose):
        print("Wrote: {}; {} planes; {:.1f} MB".format(file_out, 
                                                   len(hdu), 
                                                   os.path.getsize(file_out)/1e6))

    hdu.close()

# =============================================================================
# Make a plot if requested
# =============================================================================

    if do_plot:
        plot_backplanes(file_out, name_target, name_observer)
 
    return(0)

# =============================================================================
# End of function definition
# =============================================================================
               
# =============================================================================
# Do some tests, to validate the function
# =============================================================================

if (__name__ == '__main__'):
    
#    file_in = os.path.join(os.path.expanduser('~'), 'Data', 'ORT1', 
#                                'lor_0406731132_0x633_sci_HAZARD_test1.fit'
#                                'lor_0406731132_0x633_sci_HAZARD_test1-1.fit'
#                                'lor_0368310087_0x633_sci_HAZARD_test1.fit'
#                                )

    file_in =  '/Users/throop/Data/ORT1/porter/pwcs_ort1/K1LR_HAZ00/lor_0405175932_0x633_pwcs.fits'
    file_out = '/Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ00/lor_0405175932_0x633_pwcs_backplaned_rot.fits'
    
    name_observer = 'New Horizons'
    name_target   = 'MU69'
    frame         = '2014_MU69_SUNFLOWER_ROT'
    
    do_plot       = True

    # Start SPICE, if necessary
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')

    # Create the backplanes on disk
    
    nh_create_backplanes_fits(file_in, name_target, frame, name_observer, file_out, 
                              do_plot=True, do_verbose=True, do_clobber=True)
    
    # Plot them
    
#    plot_backplanes(file_out, name_observer = name_observer, name_target = name_target)


###
#    Now do some one-off tests. 
    
    dir = '/Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ00'
    file1 = 'lor_0405176022_0x633_pwcs_backplaned.fits'
    file2 = 'lor_0405176022_0x633_pwcs_backplaned_2.fits'
    
    hdu1 = fits.open(os.path.join(dir, file1))
    hdu2 = fits.open(os.path.join(dir, file2))
    
    r1 = hdu1['RADIUS_EQ'].data
    r2 = hdu2['RADIUS_EQ'].data
    
    plt.imshow(r1 < 100_000, alpha=0.5)
    plt.imshow(r2 < 100_000, alpha=0.5)
    plt.show()
    