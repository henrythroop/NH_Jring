#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 21:55:24 2016

@author: throop
"""

# Just a temp file to read in and examine the length of some files. Don't need to save this long-term.
# This is to check the number of frames, number of pointings, etc. relative to what Tod Lauer says.
# Concl: I agree with him.

import astropy
from   astropy.io import fits
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import hbt

hbt.figsize((60,10)) # For the long thin MVIC frames

# Load the individual frames

dir = '/Users/throop/data/NH_MVIC_Ring/'

files =  ['mpf_0299292106_0x548_sci_2.fit', # O_RingDep_A_a
          'mpf_0299793106_0x539_sci_2.fit', # O_RING_DEP_MVICFRAME_202  
          'mpf_0308706966_0x539_sci_6.fit'] # O_RING_DEP_MVICFRAME_305
          
stretch_percent = 90

hdulist = []

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 

plt.set_cmap('Greys_r')
        
for file in files:
    hdulist.append(fits.open(dir + file))

    
plt.imshow(hbt.remove_sfit(stretch(hdulist[2]['PRIMARY'].data[0,:,:]),2))

# Load the mosaics

files2 = ['ringdep_mos_v1_wcs.fits', 'mvic_d202_mos_v1.fits', 'mvic_d211_mos_v1_wcs.fits', 
                'mvic_d305_sum_mos_v1_wcs.fits']    
hdulist2 = []

for file in files2:
    hdulist2.append(fits.open(dir + file))

for hdu in hdulist2:
    print("Size = ")

plt.show()
    
#==============================================================================
# Now a few lines to read some images, and see what the DN levels are of the mosaics compared to the individual frames.
# Are the added, averaged, etc?
#==============================================================================

# Concl: I really can't tell. The noise characteristics (ie, gaussian width) change due to binning many frames.
# And the zero-bias changes (due to normalizing / stray light). So the distributions of DN look very different
# on the originals vs. mosaics, and I can't tell whether the units are ultimately the same or different.
 
hbt.figsize((8,5)) # For the long thin MVIC frames

# Make a plot for the distribution of the raw frames (pre-mosaic)

arr = hdulist[0]['PRIMARY'].data.flatten()

for i in range(np.shape(hdulist)[0]):
    arr = hdulist[i]['PRIMARY'].data
    plt.hist(arr[1].flatten(), bins=3000)
    plt.title(hdulist[i]['PRIMARY'].header['SAPNAME'] + '    ' + files[i] + ', RAW DATA')
    plt.xlim((-1,100))
    plt.show()
    print("Median = " + repr(np.median(arr)))


#for i in range(np.shape(hdulist)[0]):
#    arr = hdulist[i]['PRIMARY'].data
#    plt.hist(arr[1].flatten(), bins=10000)
#    plt.title(hdulist[i]['PRIMARY'].header['SAPNAME'] + '    ' + files[i])
#    plt.xlim((-1,100))
#    plt.show()
    
# Make a plot for distribution of mosaiced data

for i in range(np.shape(hdulist2)[0]):
    arr2 = hdulist2[i]['PRIMARY'].data   # Returns single 5000 x 700 array 

    plt.hist(arr2.flatten(), bins=2000)
    plt.title(files2[i])
    plt.xlim((-1,100))
    plt.show()
    print("Median = " + repr(np.median(arr2)))


for i in range(np.shape(hdulist2)[0]):
    arr2 = hdulist2[i]['PRIMARY'].data   # Returns single 5000 x 700 array 

    plt.hist(arr2.flatten(), bins=10000)
    plt.title(files2[i] + ", MOSAIC")
    plt.xlim((-1,10))
    plt.show()
    print("Median = " + repr(np.median(arr2)))
    
    
plt.hist(hdulist[1]['PRIMARY'].data.flatten(), bins=300)
plt.show()

plt.hist(hdulist[2]['PRIMARY'].data.flatten(), bins=300)


