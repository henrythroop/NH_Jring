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
import matplotlib.pyplot as plt
import matplotlib
import pickle # For load/save
import scipy
import scipy.misc
#import Image

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

plt.rc('image', cmap='Greys_r')
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
import scipy.misc


######

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

lun = open(filename_save, 'rb')
t = pickle.load(lun)
lun.close()
	
groupmask = (t['Desc'] == 'Jupiter ring - search for embedded moons')
t_group = t[groupmask]	

frames_med = np.array([2,3,4,5,6,7,8,9,10,11])
num_frames_med = np.size(frames_med)

frame_arr = np.zeros((num_frames_med, 1024, 1024))

do_output_png = False
 
for i in range(num_frames_med):
    f = t_group['Filename'][i] # Look up filename
    print f
    frame = hbt.get_image_nh(f,frac_clip = 1)
    frame_arr[i,:,:] = frame	
				
if (do_output_png):
    for i in range(num_frames_med):				
#        scipy.misc.toimage(frame_arr[i,:,:], cmin=0.0).save('outfile' + repr(i)+ '.jpg')
#        im = Image.fromarray(frame_arr[i,:,:])
#        im.save('outfile' + repr(i) + '.png')
        matplotlib.image.imsave('outfile' + repr(i) + '.png', hbt.remove_brightest(frame_arr[i,:,:], 0.97, symmetric=True))								

# Take the median of all of these

frame_med = np.median(frame_arr, axis=0)
frame_min = np.min(frame_arr, axis=0)
frame_max = np.max(frame_arr, axis=0)


file = t_group['Filename'][0]

image = hbt.get_image_nh(file, frac_clip=1)

plt.imshow(hbt.remove_brightest(image, 0.9, symmetric=True))
plt.title('image')
plt.show()

plt.imshow(hbt.remove_brightest(frame_med, 0.9, symmetric=True))
plt.title('median')
plt.show()

plt.imshow(hbt.remove_brightest(frame_max, 0.9, symmetric=True))
plt.title('max')
plt.show()

plt.imshow(hbt.remove_brightest(frame_min, 0.9, symmetric=True))
plt.title('min')
plt.show()

im = hbt.remove_brightest(image - frame_med, 0.9, symmetric=True)
plt.imshow(im - hbt.sfit(im,3))
plt.title('image - median, sfit')
plt.show()

im = hbt.remove_brightest(image - frame_min, 0.9, symmetric=True)
plt.imshow(im - hbt.sfit(im,3))
plt.title('image - min, sfit')
plt.show()

im = hbt.remove_brightest(image - frame_max, 0.9, symmetric=True)
plt.imshow(im - hbt.sfit(im,3))
plt.title('image - max, sfit')
plt.show()


im = hbt.remove_brightest(image - frame_arr[3,:,:], 0.9, symmetric=True)
plt.imshow(im - hbt.sfit(im,3))
plt.title('image - frame_3, sfit')
plt.show()



						
