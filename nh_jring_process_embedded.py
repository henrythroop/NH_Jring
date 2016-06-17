# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:42:18 2016

@author: throop
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 08:51:56 2016

# NH_JRING_PROCESS_EMBEDDED.PY
#
# This routine reads NH J-RING files from the 'Jupiter ring - search for embedded moons' sequence.
# It does stray light subtraction on them, and repoints, then writes them out as pngs, where they can be animated.
# Flux is preserved until the end.
#
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

######
# Now look at how to create a good median filter to subtract. 
# These routines here are unrelated to the ones above.
######

plt.rc('image', cmap='Greys_r')

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in
dir_out = '/Users/throop/data/NH_Jring/out/'

lun = open(filename_save, 'rb')
t = pickle.load(lun)
lun.close()
	
groupmask = (t['Desc'] == 'Jupiter ring - search for embedded moons')
t_group = t[groupmask]	


frames_med = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
num_frames_med = np.size(frames_med) 

frames_data = hbt.frange(0,48,49) # 0 .. 49

frame_arr      = np.zeros((num_frames_med, 1024, 1024))
frame_sfit_arr = np.zeros((num_frames_med, 1024, 1024))
frame_ffit_arr = np.zeros((num_frames_med, 1024, 1024))

power = 5

# read frames in to make a median
 
for i in range(num_frames_med):
    f = t_group['Filename'][i] # Look up filename
    print "Reading: " + f
    frame = hbt.get_image_nh(f,frac_clip = 1)
    frame_arr[i,:,:] = frame	
    
    frame_sfit_arr[i,:,:] = frame - hbt.sfit(frame, power)
    frame_ffit_arr[i,:,:] = frame - hbt.ffit(frame)

frame_sfit_med = np.median(frame_sfit_arr, axis=0) # Take median.
frame_ffit_med = np.median(frame_ffit_arr, axis=0) # Take median.

# Now read all the frames, and process them one by one

## *** Need to roll the image properly as per navigation!
 
fs = 15

for i in frames_data:
#for i in range(5):

    outfile = dir_out + 'embedded_med' + repr(int(num_frames_med)) + '_frame' + repr(int(i)) + '.png'    
    
    f = t_group['Filename'][i] # Look up filename
    print f
    print "Writing: " + f
    frame = hbt.get_image_nh(f,frac_clip = 1)
    im = hbt.remove_brightest(frame - frame_sfit_med, 0.97, symmetric=True)
    im2 = hbt.remove_brightest(im - hbt.sfit(im,8), 0.97, symmetric=True)
    
    arr_x = eval('np.' + t_group['x_pos_ring1'][i])
    arr_y = eval('np.' + t_group['y_pos_ring1'][i])
    ymean = np.mean(arr_y)
    xmean = np.mean(arr_x)
         
    dx = -int(xmean - 400)+100
    dy = -int(ymean - 3900)
    
    pad_dx = 200
    pad_dy = 200
    
    im_pad = np.zeros((1024 + pad_dx, 1024 + pad_dy))
    im_pad[pad_dx/2:-pad_dx/2, pad_dy/2:-pad_dy/2] = im2
    
    im_roll = np.roll(np.roll(im_pad, dx, axis=1), dy, axis=0)

    plt.imshow(im_roll, vmin=-20, vmax=20)
    plt.title('image ' + repr(i) + ', ' + t_group['Shortname'][i] + ', xmean=' + repr(int(xmean)) + 
               ', ymean=' + repr(int(ymean)), fontsize=fs)    
    
    matplotlib.image.imsave(outfile, im_roll)
    plt.show()
    
stop

#plt.subplot(3,2,1) # Reversed! num cols, num rows ** Good reduction. Largest az etent.
#im = hbt.remove_brightest(image - frame_sfit_med, 0.97, symmetric=True)
#plt.imshow(hbt.remove_brightest(im - hbt.sfit(im,8), 0.97, symmetric=True), vmin=-20, vmax=20)
#plt.title('image - multi_median_sfit, sfit', fontsize=fs)
##plt.show()
#
#plt.subplot(3,2,2) # OK reduction. Bit lumpy though.
#im = hbt.remove_brightest(image - frame_ffit_med, 0.97, symmetric=True)
#plt.imshow(hbt.remove_brightest(im - hbt.ffit(im), 0.97, symmetric=True), vmin=-20, vmax=20)
#plt.title('image - multi_median_ffit, ffit', fontsize=fs)
##plt.show()
#
#plt.subplot(3,2,3) # ** Good reduction. Less az extent, but maybe cleaner? 
#im = hbt.remove_brightest(image - frame_sfit_med, 0.97, symmetric=True)
#plt.imshow(hbt.remove_brightest(im - hbt.ffit(im), 0.97, symmetric=True), vmin=-20, vmax=20)
#plt.title('image - multi_median_sfit, ffit', fontsize=fs)
##plt.show()
#
#plt.subplot(3,2,4) # ** Not good. Too lumpy.
#im = hbt.remove_brightest(image - frame_ffit_med, 0.97, symmetric=True)
#plt.imshow(hbt.remove_brightest(im - hbt.sfit(im,8), 0.97, symmetric=True), vmin=-20, vmax=20)
#plt.title('image - multi_median_ffit, sfit', fontsize=fs)
#
## Now do two more, denoised. For these, use a median filter but only in the vertical direction (10,1).
## Otherwise the median will lose its vertical zebra stripes. Those are a source of noise that 
## we actually want to keep in the median, so they are removed when we subtract them.
#
#plt.subplot(3,2,5) # *** I think this is the best overall reduction
#frame_sfit_med_denoise = scipy.ndimage.median_filter(frame_sfit_med,size=(10,1))
#im = hbt.remove_brightest(image - frame_sfit_med_denoise, 0.97, symmetric=True)
#im2 = hbt.remove_brightest(im - hbt.ffit(im), 0.97, symmetric=True)
#plt.imshow(im2, vmin=-20, vmax=20)
#
#plt.title('image - multi_median_sfit_denoise, ffit', fontsize=fs)
#
#plt.subplot(3,2,6) # *** too lumpy
#frame_ffit_med_denoise = scipy.ndimage.median_filter(frame_ffit_med,size=(10,1))
#im = hbt.remove_brightest(image - frame_ffit_med_denoise, 0.97, symmetric=True)
#plt.imshow(hbt.remove_brightest(im - hbt.sfit(im,8), 0.97, symmetric=True), vmin=-20, vmax=20)
#plt.title('image - multi_median_ffit_denoise, sfit', fontsize=fs)
#plt.show()



