# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 09:37:07 2016

NH_JRING_CREATE_MEDIANS.PY

This routine creates median files, used for stray light subtraction.

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
import os

#import Image

######
# Various routines to create the straylight median files
######

def nh_get_straylight_median(index_group, index_files, do_fft=False, do_sfit=True, power1=5, power2=5):
    
    """
    This is a user-level routine to get a straylight median image.
    If the image exists, it will load it. If it doesn't exist, it will create it.
    This is the one to call under most circumstances.
    """
    
    dir_straylight = '/Users/throop/data/NH_Jring/out/'
    
    file_base = nh_create_straylight_median_filename(index_group, index_files, 
                                                      do_fft=do_fft, do_sfit=do_sfit, power1=power1, power2=power2)
    
    file_pkl = dir_straylight + file_base + '.pkl'
    
    # If file exists, load it and return
    
    if os.path.exists(file_pkl):
#           print 'Reading file: ' + file_pkl
        lun = open(file_pkl, 'rb')
        arr = pickle.load(lun)
        lun.close() 
        return arr

    # If file doesn't exist, create it, save it, and return
        
    else:
        (arr, file) = nh_create_straylight_median(index_group, index_files, 
                                                      do_fft=do_fft, do_sfit=do_sfit, power1=power1, power2=power2)

        lun = open(file_pkl, 'wb')
        pickle.dump(arr, lun)
        lun.close()
        print 'Wrote file: ' + file_pkl
        
        return arr
        
def nh_create_straylight_median(index_group, index_files, do_fft=False, do_sfit=True, power1=5, power2=5):
    
    """
    This takes a set of related observations, and creates a median image of all of them,
    in a form useful for straylight removal.
    For now this routine is hard-coded to NH J-ring, but it could be generalized later.
    """

#     o group:   What set of observations, grouped by 'Desc' field. Integer.
#     o files:   Int array list of the files 
#     o do_fft:  Flag. For the final step, do we use an FFT or sfit?
#     o do_sfit: Flag.
#     o power1:  Exponent for sfit, to be applied to each frame (step 1)
#     o power2:  Exponent for sfit, to be applied to each frame (step 2)
#   
#     This routine returns the array itself, and a recommended base filename. It does not write it to disk.
 
    file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
    
    lun = open(file_pickle, 'rb')
    t = pickle.load(lun)
    lun.close()

    # Process the group names. Some of this is duplicated logic -- depends on how we want to use it.

    groups = astropy.table.unique(t, keys=(['Desc']))['Desc']
    
#    groupname = 'Jupiter ring - search for embedded moons'
#    groupnum = np.where(groupname == groups)[0][0]
        
    groupmask = (t['Desc'] == groups[index_group])
    t_group = t[groupmask]	
    
    # Create the output arrays
    # For now I am assuming 1X1... I'll have to rewrite this if I use 4X4's very much.
   
    num_files = np.size(index_files)
    
    frame_arr      = np.zeros((num_files, 1024, 1024))
    frame_sfit_arr = np.zeros((num_files, 1024, 1024))
    frame_ffit_arr = np.zeros((num_files, 1024, 1024))
    
    # Read frames in to make a median
     
    for i,n in enumerate(index_files):
        file = t_group['Filename'][n] # Look up filename
        print "Reading: " + file
        frame = hbt.read_lorri(file,frac_clip = 1)
        if (np.shape(frame)[0] == 256):
            
    # Resize the image to 1024x1024, if it is a 4x4 frame. 
    # scipy.misc.imresize should do this, and has various interpolation schemes, but truncates to integer (?!)
    # because it is designed for images. So instead, use scipy.ndimage.zoom, which uses cubic spline.
        
            frame2 = scipy.ndimage.zoom(frame, 4)
            frame = frame2
        frame_arr[i,:,:] = frame	
        
        frame_sfit_arr[i,:,:] = frame - hbt.sfit(frame, power1)
        frame_ffit_arr[i,:,:] = frame - hbt.ffit(frame)
    
    frame_sfit_med = np.median(frame_sfit_arr, axis=0) # Take median using sfit images
    frame_ffit_med = np.median(frame_ffit_arr, axis=0) # Take median using fft  images
    
    # Do a very small removal of hot pixels. Between 0.99 and 0.999 seems to be fine. Make as small 
    # as possible, but something is often necessary
    
#    frame_sfit_med = hbt.remove_brightest(frame_sfit_med, 0.999)
#    frame_ffit_med = hbt.remove_brightest(frame_sfit_med, 0.999)
#    
    file_base = nh_create_straylight_median_filename(index_group, index_files, do_fft=do_fft, do_sfit=do_sfit, 
                                                     power1=power1, power2=power2)
                                                     
    if (do_fft):   
        return (frame_ffit_med, file_base)
    
    if (do_sfit):
        return (frame_sfit_med, file_base)
        
def nh_create_straylight_median_filename(index_group, index_files, do_fft=False, do_sfit=True, power1=5, power2=5):
    """
    Just create the base filename for a stray likght median file.
    """
    # Usually we will add extensions onto this:
    #  .pkl -- save file
    #  .png -- file to be edited in photoshop
    
    if (do_fft):    #fig, ax = plt.subplots(2, 1, figsize=(20, 15 ))
        file_base = 'straylight_median_g{:.0f}_n{:.0f}..{:.0f}_fft'.format(index_group, np.min(index_files), np.max(index_files))
        return file_base
    
    # NB: The print format string :.0f is useful for taking a float *or* int and formatting as an int.
    #     The 'd' string is for a decimal (int) but chokes on a float.
    
    if (do_sfit):
        file_base = 'straylight_median_g{:.0f}_n{:.0f}..{:.0f}_sfit{:.0f},{:.0f}'.format(index_group, 
                                                             np.min(index_files), np.max(index_files), power1, power2)
        return file_base
        
        
def nh_create_straylight_medians():

    """
    This is a wrapper routine which calls the main function. It is just for testing. 
    """
    
    plt.rc('image', cmap='Greys_r')
     
    segment = 4  # 1, 2, or 3
    
    dir_out = '/Users/throop/data/NH_Jring/out/' 
    
    fs = 15 # Set the font size
    
    if (segment == 1): # 50 frames, upright
        frames_med = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
        frames_med = np.array([1,2,3,4,5])
        index_files = frames_med
        index_group = 8
        
    #    num_frames_med = np.size(frames_med)     
    
        frames_data = hbt.frange(0,48,49) # 0 .. 49
    
    if (segment == 2): # 5 frames, tilted sideways
        frames_med = np.array([49,51,52,53,53])
        num_frames_med = np.size(frames_med) 
        frames_data = frames_med
    
    if (segment == 3): # 5 frames, tilted sideways
        frames_med = hbt.frange(54,111,58)
        num_frames_med = np.size(frames_med) 
        frames_data = frames_med

    # Get the straylight frame

    im = nh_create_straylight_median(index_group, index_files, do_fft=False, do_sfit=True, power1=5, power2=5)

    # Save the straylight frame itself
    
    outfile = dir_out + 'embedded_seg' + repr(int(segment)) + '_med' + repr(int(num_frames_med)) + '_sfit.png'
    matplotlib.image.imsave(outfile, frame_sfit_med)
    print 'Wrote: ' + outfile
    
    outfile = dir_out + 'embedded_seg' + repr(int(segment)) + '_med' + repr(int(num_frames_med)) + '_ffit.png'
    matplotlib.image.imsave(outfile, frame_ffit_med)
    print 'Wrote: ' + outfile
    
    file_out = 'med_g9_n0-49_sfit8_571_ps.pkl'
    
#for i,n in enumerate(frames_data):
#for i in range(5):

    outfile = dir_out + 'embedded_seg' + repr(int(segment)) + '_med' + repr(int(num_frames_med)) + \
        '_frame' + repr(int(i)) + '.png'    
    
    file = t_group['Filename'][n] # Look up filename
    print file
    print "Writing: " + outfile
    frame = hbt.read_lorri(file,frac_clip = 1)
    im  = hbt.remove_brightest(frame - frame_sfit_med, 0.97, symmetric=True)
    im2 = hbt.remove_brightest(im    - hbt.sfit(im,power2), 0.97, symmetric=True)
    
    arr_x = eval('np.' + t_group['x_pos_ring2'][n]) # Ring2 = outer ring
    arr_y = eval('np.' + t_group['y_pos_ring2'][n])
    ymean = np.mean(arr_y)
    xmean = np.mean(arr_x)
    
    ymin  = np.min(arr_y) # top of screen -- for upward pointing ansae
    ymax  = np.max(arr_y) # Bottom of screen -- for downward pointing ansae
    xmin  = np.min(arr_x) # Left edge
    xmax  = np.max(arr_x) # Right edge
    
    if (segment == 1): #     
        roll_x = -int(xmin - 100)    # Roll to __ from left edge
        roll_y = -int(ymin - 200)    # Roll to __ from top  edge
    
    if (segment == 2):
        roll_x = -int(xmean - 400)+500
        roll_y = -int(ymean - 3900)-1050

    if (segment == 3):     
        roll_x = -int(xmean - 400)+100
        roll_y = -int(ymean - 2600)    # -3200 is toward top. -3400 is toward bottom 


    # Now read all the frames, and process them one by one.
    # Output as images which can be animated.
    
    ## *** Need to roll the image properly as per navigation!

         
    pad_dx = 400 # Amount to add in x dir, total
    pad_dy = 400
    
    im_pad = np.zeros((1024 + pad_dx, 1024 + pad_dy))-20
    im_pad[pad_dx/2:-pad_dx/2, pad_dy/2:-pad_dy/2] = im2
    
    im_roll = np.roll(np.roll(im_pad, roll_x, axis=1), roll_y, axis=0)

    ax = plt.imshow(im_roll, vmin=-20, vmax=20)
    plt.title(repr(i) + ', ' + t_group['Shortname'][n])
    plt.axis('off')
    
#+ ', xmean=' + repr(int(xmean)) + 
#               ', ymean=' + repr(int(ymean)), fontsize=fs)    
    
    # Save the image
    
    matplotlib.image.imsave(outfile, im_roll)
    plt.show()

    plt.subplot(1,2,1) # Reversed! num cols, num rows ** Good reduction. Largest az etent.
    plt.imshow(frame_sfit_med)
    plt.title('frame_sfit_med, n=' + repr(num_frames_med))
    #plt.show()
    
    plt.subplot(1,2,2) # Reversed! num cols, num rows ** Good reduction. Largest az etent.
    plt.imshow(frame_ffit_med)
    plt.title('frame_ffit_med, n=' + repr(num_frames_med))
    plt.show()    
    
#stop




#####

#def f():
#    nh_get_straylight_median(8,[50,51,52,52], do_sfit=True, power1=5, power2=5) # long movie, tilted portion
#    nh_get_straylight_median(8,hbt.frange(0,49,50), do_sfit=True, power1=5, power2=5) # long movie, first 50 frames
#nh_get_straylight_median(5,hbt.frange(1,6,6), do_sfit=True, power1=5, power2=5) # best portrait, ansa

s = nh_get_straylight_median(6,hbt.frange(229,232), do_sfit=True, power1=5, power2=5) # 4x4 gossamer 
a = hbt.read_lorri(35117774)
a = hbt.read_lorri(35120054)
a = a - hbt.sfit(a,5)
a = hbt.remove_brightest(a,0.99,symmetric=True)  # Careful here -- don't crop the ring!
(s_norm,coeffs) = hbt.normalize_images(s,a)
plt.imshow(a - s_norm)

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



