# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 08:51:56 2016

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

def remove_stray():
    
    return 0
    
dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
file1 = 'lor_0034602723_0x630_sci_1.fit' # 'embedded moons', position 1. Pointed down.
file2 = 'lor_0034603323_0x630_sci_1.fit' # 'embedded moons', position 2. Pointed down.
file3 = 'lor_0034676524_0x630_sci_1.fit' # 'phase curve', all fixed position. Pointed up.
file4 = 'lor_0034676944_0x630_sci_1.fit' # 'phase curve', all fixed position. Pointed down.

# Look up the sequence and ET for these images!

hdulist = fits.open(dir_images + file1)
    
print hdulist[0].header['REQDESC'] 
print hdulist[0].header['SPCUTCAL']
hdulist.close()
    
a1 = hbt.get_image_nh(dir_images + file1, bg_method='raw', frac_clip = 1.)
a2 = hbt.get_image_nh(dir_images + file2, bg_method='raw', frac_clip = 1.)
a3 = hbt.get_image_nh(dir_images + file3, bg_method='raw', frac_clip = 1.)
a4 = hbt.get_image_nh(dir_images + file4, bg_method='raw', frac_clip = 1.)

hbt.set_plot_defaults()

#plt.rc('image', cmap='Greys')
plt.rcParams['figure.figsize'] = (10, 10)

plt.imshow(a1)
plt.show()

plt.imshow(a2)
plt.show()

# Reproduce what Hal and Andy did

d = a1 - 0.756 * a2
d2 = -(hbt.remove_brightest(-hbt.remove_brightest(d, 0.95), 0.95)) # Remove both the brightest... and darkest!
plt.imshow(d2)
plt.show()

# 

plt.imshow(hbt.remove_brightest(a3, 0.9))
plt.show()

diff = -(hbt.remove_brightest(-hbt.remove_brightest(a4 - 0.5*a2, 0.97), 0.97))
diff_s = diff - hbt.sfit(diff, degree=5)
plt.imshow(diff_s)
plt.ylim((1023, 0))
plt.show()
