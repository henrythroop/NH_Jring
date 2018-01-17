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

from nh_create_backplanes_fits import nh_create_backplanes_fits

def nh_make_backplanes_ort1():
    
    """
    Process all of the MU69 ORT1 files. 
    
    Takes Simon's WCS'd FITS files as inputs, and creates backplaned FITS as output.
    
    Call this function in order to generate the backplanes from Simon's WCS files.
    """
    
# =============================================================================
#     Get a proper list of all the input files
# =============================================================================
    
    dir_data_ort1 = '/Users/throop/Data/ORT1'
    dir_in  = os.path.join(dir_data_ort1, 'porter', 'pwcs_ort1')
    dir_out = os.path.join(dir_data_ort1, 'throop', 'backplaned')
        
    do_plot    = False
    do_clobber = False
    do_reverse = True
    do_digit_filter = False
    
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
            digit = base[12]  # This is the digit that changes the most in the LORRI files, so it's a good choice.
            if (digit in digit_filter):
                files_filtered.append(file)
        print("Filtered on '{}': {} files → {}".format(digit_filter, len(files), len(files_filtered)))
        
        files = files_filtered            

# =============================================================================
# Now that we have a list of files, go thru and process each one.
# =============================================================================
        
    for i,file_in in enumerate(files):
        print("{}/{}".format(i,len(files))) 
        file_out = file_in.replace(dir_in, dir_out)
        file_out = file_out.replace('_pwcs.fit', '_pwcs_backplaned.fit') # Works for both .fit and .fits
    
        try:
            nh_create_backplanes_fits(file_in, 
                                      file_out, 
                                      do_plot=False, 
                                      do_clobber=do_clobber,
                                      verbose=True)
        except FileExistsError:
            print('File exists -- skipping. {}'.format(file_out))
    
        if (do_plot):
            nh_plot_backplanes(file_out)
 
# =============================================================================
# End of function
# =============================================================================
    
# =============================================================================
# Run the function if requested
# =============================================================================
           
if (__name__ == '__main__'):
    nh_make_backplanes_ort1()
            