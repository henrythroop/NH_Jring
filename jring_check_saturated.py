# -*- coding: utf-8 -*-
"""

This is a one-off program to check saturation levels on a couple of NH LORRI frames.

Created on Thu Jul  7 12:50:49 2016

@author: throop
"""

import math      
import astropy
from   astropy.io import fits
import numpy as np
import cspice
import wcsaxes
import hbt
from   astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib
import pickle # For load/save
import scipy
import scipy.misc

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'

# We read level-2 images, which return 'photometrically corrected DN', as per ICD @ 44
# These correct for flat-field and bias, but not stray light etc. See ICD.

file1 = 'lor_0035123700_0x633_sci_1.fit'  # 4x4 LORRI image. Typical DN level 3600 [saturated, makes sense]
file2 = 'lor_0035282100_0x633_sci_1.fit'  # 4x4 LORRI image. Typical DN level -1.5 [saturated, but not sure why neg]

file = dir_images + '/' + file1

hdulist = fits.open(file)

hdu_image = hdulist['PRIMARY'] # Options are 'PRIMARY', 'LORRI Error image', 'LORRI Quality flag image'. 
                               # Get list with hdulist.info(), which is different than hdulist.info, which is also valid.

hdu_error = hdulist['LORRI Error image']
hdu_quality = hdulist['LORRI Quality flag image']

# New Horizons ICD (in papers.app) @ 45: 16 = "Saturated pixel in raw data (A/D value of 4095)"

image = hdu_image.data
quality = hdu_quality.data
error = hdu_error.data

hdulist.close()

image_proc = image - hbt.sfit(image,5)
image_proc = hbt.remove_brightest(image_proc,0.9,symmetric=True)
plt.imshow(image_proc)
plt.title(file)
plt.show()

# Now print summary of the images.
#
# Looks like the 4x4 image ('23700') has DN level of 3500 with 10 sec
# 1x1 image ('0624') has DN level of 123 with 3 sec.
# So, we would expect 4x4 to have 16 x 10 / 3 ~ 50x the DN level of the 1x1. And indeed that is true.
# So images are bascially saturated as expected. Ugh. 
# The '2100' image has a negative DN value. Not sure why. But it is also a long 4x4 - intended for depth -- and
# just got swamped by stray light. Not by a super amount... but just barely. Would have been just fine as 1x1's, although
# maybe would have had a small amount of smear.
# This agrees w/ e-mail from Hal 3-Jul-2016. But I just wanted to load images and verify.

print "file = " + file
print "mean quality = " + repr(np.mean(quality))
print "mean DN = " + repr(np.mean(image))
print "mean error = " + repr(np.mean(error))
    
print
    
file3 = 'lor_0035120624_0x630_sci_1.fit'
image_1x1 = hbt.get_image_nh(dir_images + '/' + file3)
print "file = " + file3
print "mean DN = " + repr(np.mean(image_1x1))
