#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 23:42:14 2018

@author: throop
"""

import math      
import astropy
from   astropy.io import fits
import numpy as np
import spiceypy as sp
from   astropy.visualization import wcsaxes
import hbt
from   astropy.wcs import WCS
import os
import matplotlib.pyplot as plt
import glob

from create_backplanes_fits import create_backplanes_fits
from plot_backplanes        import plot_backplanes

def nh_ort_make_backplanes(frame = '2014_MU69_SUNFLOWER_ROT', digit_filter=None):
    
    """
    Process all of the MU69 ORT files. 
    
    Takes Simon's WCS'd FITS files as inputs, and creates backplaned FITS as output.
    
    Call this function in order to generate the backplanes from Simon's WCS files.
    """

# =============================================================================
# Initialize
# =============================================================================

    do_plot    = True
    do_clobber = True
    
    name_target   = 'MU69'
    name_observer = 'New Horizons'
    # frame         = '2014_MU69_SUNFLOWER_ROT'  # Change this to tuna can if needed, I think??
    # frame         = '2014_MU69_TUNACAN_ROT'
    # frame         = '2014_MU69_ORT4_1'  # Change this to tuna can if needed, I think??
    
# =============================================================================
#     Get a proper list of all the input files
# =============================================================================
    
    do_ORT1 = False
    do_ORT3 = False
    do_ORT2 = False
    do_ORT4 = False
    do_ACTUAL = True  # Run this on actual OpNav data!
    
    do_force = True
    
#    dir_data_ort = '/Users/throop/Data/ORT1'
#    dir_in  = os.path.join(dir_data_ort, 'porter', 'pwcs_ort1')
#    dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')

    if do_ORT2:
        dir_data_ort = '/Users/throop/Data/ORT2'
        dir_in  = os.path.join(dir_data_ort, 'porter', 'pwcs_ort2')
        dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')
        files = glob.glob(os.path.join(dir_in, '*','*_ort2.fit'))

    if do_ORT3:    
        dir_data_ort = '/Users/throop/Data/ORT3'
        dir_in  = os.path.join(dir_data_ort, 'buie') # Using Buie backplanes, not Simon's.
        dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')
        files = glob.glob(os.path.join(dir_in, '*','*_ort3.fit'))
        frame         = '2014_MU69_TUNACAN_ROT'
        
    if do_ORT4:
        dir_data_ort = '/Users/throop/Data/ORT4'
        dir_in  = os.path.join(dir_data_ort, 'porter', 'pwcs_ort4')
        dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')
        files = glob.glob(os.path.join(dir_in, '*','*_pwcs.fits'))

    if do_ACTUAL:
        dir_data_ort = '/Users/throop/Data/MU69_Approach'
        dir_in  = os.path.join(dir_data_ort, 'porter')
        dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')
        files = glob.glob(os.path.join(dir_in,'*', '*_pwcs.fits'))  # OpNav field data (older)
        
        files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018267', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018284', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018284', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018284', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018284', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018287', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018298', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018301', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018304', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018311', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018314', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018315', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018316', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018317', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018325', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018326', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018327', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_OpNav_L4_2018328', '*_pwcs2.fits')) # OpNav data
        # files = glob.glob(os.path.join(dir_in, 'KALR_MU69_Hazard_L4_2018325', '*_pwcs2.fits')) # OpNav data

        files = glob.glob(os.path.join(dir_in, '*', '*_pwcs2.fits')) # OpNav data
        
        do_force = False
        do_clobber = False

# =============================================================================
# Check what FRAME we are using, and change output directory appropriately
# =============================================================================

    if 'TUNACAN' in frame.upper():
        dir_out = dir_out.replace('backplaned', 'backplaned_tunacan')
        
# =============================================================================
#     Filter files if needed
# =============================================================================
    
    # If desired, do a 'digit filter.' This filters the files down into a smaller number.
    # This is useful to do processing in parallel. Python global interpreter lock means
    # that only one CPU at a time can be used. To get around this, filter the files down,
    # and put each filter in its own Spyder tab.
    
    if digit_filter:
        
    # do_digit_filter = False

    # digit_filter = '12'
    # digit_filter = '34'
    # digit_filter = '56'
    # digit_filter = '78'
    # digit_filter = '90'
    
    # if (do_digit_filter):
        files_filtered = []
        for file in files:
            base = os.path.basename(file)
            digit = base[12]  # Match the penultimate digit in LORRI filename (e.g., lor_0405348852_pwcs ← matches '5')
                              # This is the digit that changes the most in the LORRI files, so it's a good choice.
            if (digit in digit_filter):
                files_filtered.append(file)
        print("Filtered on '{}': {} files → {}".format(digit_filter, len(files), len(files_filtered)))
        
        files = files_filtered            

# =============================================================================
# Start SPICE, if necessary
# =============================================================================
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
        
# =============================================================================
# Loop and create each backplane
# =============================================================================
        
    for i,file_in in enumerate(files):
        print("{}/{}".format(i,len(files))) 
        file_out = file_in.replace(dir_in, dir_out)
        file_out = file_out.replace('_pwcs.fit', '_pwcs_backplaned.fit') # Works for both .fit and .fits
        file_out = file_out.replace('_pwcs2.fit', '_pwcs2_backplaned.fit') # Works for both .fit and .fits
    
        # Call the backplane function. Depending on settings, this will automatically run if a newer input file is 
        # received, and thus we need to regenerate the output backplane.
        
        try:
            create_backplanes_fits(file_in, 
                                      name_target,
                                      frame,
                                      name_observer,
                                      file_out,
                                      do_plot=False, 
                                      do_clobber=do_clobber,
                                      do_verbose=True)
            if (do_plot):
                plot_backplanes(file_out, name_observer = name_observer, name_target = name_target)
     
        except FileExistsError:
            print('File exists -- skipping. {}'.format(os.path.basename(file_out)))
    
       
# =============================================================================
# End of function
# =============================================================================

def test():

    dir = '/Users/throop/Data/ORT3/throop/backplaned'
    file = 'lor_0406991502_0x633_wcs_HAZARD_ort3_backplaned.fit'
    file = pwcs_ort4/K1LR_OPNAV27B/lor_0407015627_0x633_pwcs.fits
    plane = nh_ort_make_backplanes(os.path.join(dir,file))
    
# =============================================================================
# Run the function if requested
# When run, this program regenerates all of the backplanes    
# =============================================================================
        
if (__name__ == '__main__'):
    file_tm = 'kernels_kem_prime.tm'
    sp.unload(file_tm)
    sp.furnsh(file_tm)

 # NB: This would be an ideal candidate for multi-processing
        
    # Set paramters here
    
    do_tuna = False
    digit_filter = None
    
    # digit_filter = '12'
    # digit_filter = '34'
    # digit_filter = '56'
    # digit_filter = '78'
    # digit_filter = '90'
    
    # Run code here

    if do_tuna:
        print("WARNING: USING TUNACAN FRAME!")
        nh_ort_make_backplanes(frame='2014_MU69_TUNACAN_ROT', digit_filter=digit_filter)
        print("WARNING: USING TUNACAN FRAME!")
        
    else:
        nh_ort_make_backplanes(digit_filter=digit_filter)
        
        