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

def nh_make_backplanes_ort1():
    
    """
    Process all of the MU69 ORT1 files. 
    
    Takes Simon's WCS'd FITS files as inputs, and creates backplaned FITS as output.
    
    Call this function in order to generate the backplanes from Simon's WCS files.
    """

# =============================================================================
# Initialize
# =============================================================================

    do_plot    = False
    do_clobber = False
    do_digit_filter = False
    
    name_target   = 'MU69'
    name_observer = 'New Horizons'
    frame         = '2014_MU69_SUNFLOWER_ROT'
    
# =============================================================================
#     Get a proper list of all the input files
# =============================================================================
    
    do_ORT1 = False
    do_ORT2 = True
    
#    dir_data_ort = '/Users/throop/Data/ORT1'
#    dir_in  = os.path.join(dir_data_ort, 'porter', 'pwcs_ort1')
#    dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')

    dir_data_ort = '/Users/throop/Data/ORT2'
    dir_in  = os.path.join(dir_data_ort, 'porter', 'pwcs_ort2')
    dir_out = os.path.join(dir_data_ort, 'throop', 'backplaned')
    
    files = glob.glob(os.path.join(dir_in, '*','*_pwcs.fits'))
    
# =============================================================================
#     Filter files if needed
# =============================================================================
    
    # If desired, do a 'digit filter.' This filters the files down into a smaller number.
    # This is useful to do processing in parallel. Python global interpreter lock means
    # that only one CPU at a time can be used. To get around this, filter the files down,
    # and put each filter in its own Spyder tab.
    
#    digit_filter = '12'
#    digit_filter = '34'
#    digit_filter = '56'
#    digit_filter = '78'
    digit_filter = '90'
    
    if (do_digit_filter):
        files_filtered = []
        for file in files:
            base = os.path.basename(file)
            digit = base[12]  # Match the penultimate digit in LORRI filename (e.g., lor_0405348852_pwcs  ← matches '5')
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
            print('File exists -- skipping. {}'.format(file_out))
    
       
# =============================================================================
# End of function
# =============================================================================
    
# =============================================================================
# Run the function if requested
# =============================================================================
           
if (__name__ == '__main__'):
    nh_make_backplanes_ort1()
            