# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 22:45:22 2016

@author: throop
"""

# Create all backplanes. This is a one-off script that is run (in theory) only one time ever.
# It creates all of the navigation backplanes based on WCS keywords found in the image headers.
# Any offset (from stellar navigation, or user-determined coeffs) is *not* added, unless the image
# has had its WCS coords updated (e.g. *_opnav.fit)
# 
# This routine takes ~ 1 minute per file.
#
# I should investigate if I can multi-thread this. It seems like it could be sped up very easily.

import math      
import astropy
from   astropy.io import fits
import numpy as np
import spiceypy as sp
from   astropy.visualization import wcsaxes
import hbt
from   astropy.wcs import WCS
import pickle # For load/save
import os

import pdb
import glob
import math

DO_OVERWRITE = False

# Start up SPICE

#file_tm = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'  # SPICE metakernel
file_tm = 'kernels_nh_jupiter.tm'
sp.furnsh(file_tm)
        
# Get the full list of files

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

dir_out = '/Users/throop/data/NH_Jring/out/'

#file_list = glob.glob(dir_images + '/*fit')
#file_list = glob.glob(dir_images + '/*0034620123*1*opnav.fit')  # Navigate the opnav'd files
file_list = glob.glob(dir_images + '/*_sci_?_opnav.fit')  # Navigate the opnav'd files

files = np.array(file_list)

DO_TEST = False

if DO_TEST:
    files = files[0:3]
    
for i,file in enumerate(files):

    file_short = file.split('/')[-1]
    file_out = dir_out + file_short    
    file_out = file_out.replace('.fit', '_planes.pkl')
    
    print("{}/{}: Generating backplane for {}".format(i, np.size(files), file_short))
    
    if os.path.isfile(file_out) and not(DO_OVERWRITE):
        print("  ** File already exists! Skipping. {}".format(file_short))
    else:    
        plane = hbt.create_backplane(file)

    #    print "RA, Dec mean = " + repr(np.mean(plane['RA']) * hbt.r2d) + ', ' + repr(np.mean(plane['Dec']) * hbt.r2d)
        
    #    print "Ang_Metis meaan = " + repr(np.mean(plane['Ang_Metis']) * hbt.r2d) + ' deg'
    #    print "Ang_Adrastea mean = " + repr(np.mean(plane['Ang_Adrastea']) * hbt.r2d) + ' deg'
    #    print "Ang_Thebe mean = " + repr(np.mean(plane['Ang_Thebe']) * hbt.r2d) + ' deg'
        
        # Write one variable to a file        
    
        lun = open(file_out, 'wb')
        pickle.dump(plane, lun)
        lun.close()
    
        print("Wrote: " + file_out)
        print
