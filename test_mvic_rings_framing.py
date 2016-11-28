#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 21:55:24 2016

@author: throop
"""

# Just a temp file to read in and examine the length of some files. Don't need to save this long-term.

import astropy
from   astropy.io import fits
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np

dir = '/Users/throop/data/NH_MVIC_Ring/'
files =  ['mpf_0299292106_0x548_sci_2.fit', 'mpf_0299793106_0x539_sci_2.fit', 'mpf_0308706966_0x539_sci_6.fit']
stretch_percent = 90
hdulist = []

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 
plt.set_cmap('Greys_r')
        
for file in files:
    hdulist.append(fits.open(dir + file))

plt.imshow(hbt.remove_sfit(stretch(hdulist[2]['PRIMARY'].data[0,:,:]),2))

files2 = ['ringdep_mos_v1_wcs.fits', 'mvic_d202_mos_v1.fits', 'mvic_d211_mos_v1_wcs.fits', 
                'mvic_d305_sum_mos_v1_wcs.fits']    
hdulist2 = []

for file in files2:
    hdulist2.append(fits.open(dir + file))

for hdu in hdulist2:
    print("Size = 