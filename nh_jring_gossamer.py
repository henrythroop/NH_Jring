#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:39:08 2018

@author: throop
"""

# General python imports

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
from   subprocess import call
import warnings
import pdb
import os.path
import os
import subprocess

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astroquery.vo_conesearch import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

import scipy

# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
from tkinter import ttk
from tkinter import messagebox
tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

from   astropy.stats import sigma_clip

# HBT imports

import hbt

#%%%
plt.set_cmap('plasma')

rj_km = 71492

# Now list all of the Gossamer observations. I have determined these groupings manually, based on timings.
# In general it looks like usually these are stacks of four frames at each pointing ('footprint')

index_group = 6

index_image_list = [
#                    hbt.frange(46,48), # 1X1. Main ring. Maybe gossamer too but not interested.
#                    hbt.frange(49,53), # 1X1. Main ring. Maybe gossamer too but not interested.
#                    np.array([54]),    # 1X1. Main ring. Maybe gossamer too but not interested.
                    hbt.frange(59,62),  # Right ansa. I think Gossamer is there but hidden -- need to optimize.
                    hbt.frange(63,66),  # Left ansa. Furthest out, no ring. Use for subtraction.
                    hbt.frange(67,70),  # Left ansa. Far out, no ring. Closer than above.
                    hbt.frange(71,74),  # Left ansa. No ring. Closer than above.
                    hbt.frange(75,78),  # Left ansa. Gossamer. Closer than above.
#                    hbt.frange(79,88),  # Left ansa. Main ring. Prob gossamer too but not interested. 1X1. Svrl ptgs.
#                    hbt.frange(89,94),  # Io! Closest of the sequence, and I will ignore. 1X1.
                    hbt.frange(95,98),  # Right ansa. Gossamer limb.
                    hbt.frange(99,102), # Left ansa. Similar geometry as 75-78.
                    hbt.frange(112,115),# Right ansa. Gossamer limb. Similar geometry as 95-98 (offset pointing)
                    hbt.frange(116,119),# Left ansa. Similar geometry as 75-78. Sat visible??
                    hbt.frange(120,123),  # Left ansa. Closer than above.
                    hbt.frange(124,127),  # Left ansa. Closer than above.
                    hbt.frange(128,131),  # Left ansa + main ring ansa. Closer than above.
                    hbt.frange(157,160),
                    hbt.frange(173,176),
                    hbt.frange(177,180),
                    hbt.frange(181,184),
                    hbt.frange(185,188),
                    hbt.frange(207,210),
                    hbt.frange(211,214),
                    hbt.frange(225,227),
                    np.array([228]),
                    hbt.frange(229,232),
                    hbt.frange(233,238),
                    
                    ]
                    


do_test = False

if do_test:
    
    index_image_list = [
                    hbt.frange(99,102), # Left ansa. Similar geometry as 75-78.
                    hbt.frange(112,115),# Right ansa. Gossamer limb. Similar geometry as 95-98 (offset pointing)
                    hbt.frange(116,119),# Left ansa. Similar geometry as 75-78. Sat visible??
                    hbt.frange(120,123),  # Left ansa. Closer than above.
                    hbt.frange(124,127),  # Left ansa. Closer than above.
                    hbt.frange(128,131)]  # Left ansa + main ring ansa. Closer than above.



num_footprints = len(index_image_list)
  
#index_image = hbt.frange(71,74)
#index_image = hbt.frange(75,78)
#
#index_image = hbt.frange(95,99)
#
#index_image = hbt.frange(112,115)
#index_image = hbt.frange(116,119)
#index_image = hbt.frange(120,123)
#index_image = hbt.frange(124,127)
#index_image = hbt.frange(128,131)

#index_image = hbt.frange(100,102)
#index_image = hbt.frange(49,54)

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.


filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters 

dir_out    = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
dir_backplanes = '/Users/throop/data/NH_Jring/out/'
dir_lauer  = dir_out.replace('/out/', '/lauer/')

lun = open(dir_out + filename_save, 'rb')
t_all = pickle.load(lun)
lun.close()

abcorr = 'LT+S'

groups = astropy.table.unique(t_all, keys=(['Desc']))['Desc']

groupmask = t_all['Desc'] == groups[index_group]

t_group = t_all[groupmask]  # 
        
# =============================================================================
# Make a plot showing the center position of a lot of these frames, on plot per footprint
# =============================================================================

file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel

sp.furnsh(file_tm)

hbt.figsize((9,9))

index_group = 6
index_images = np.array(
              list(hbt.frange(10,15)) +
              list(hbt.frange(46,54)) +
              list(hbt.frange(59,135)) )

imagenum = 0

#index_images = hbt.frange(49,54)

#index_images = np.array([46,47,48])

# Set up a a bunch of arrays (really lists) to save parameters from each image. We then output these.

ra_arr      = []
dec_arr     = []
dist_proj_rj_arr = []      # Projected distance in RJ, for each image
dist_proj_rj_foot_arr = [] # Projected distance in RJ, for each footprint (ie, set of four images)
range_rj_arr = []
et_arr      = []
dx_arr      = []
dy_arr      = []
dz_arr      = []
corner_arr  = []
format_arr    = []
exptime_arr = []
name_limb_arr = []
name_limb_foot_arr = []
file_arr    = []
phase_arr   = []
utc_arr     = []
index_hbt_arr = []
index_footprint_arr = []
im_process_arr = []   # Array of final processed images
im_lauer_arr = []
im_sum_s_arr = []
    
for index_footprint,index_images in enumerate(index_image_list):  # Loop over *all* the images
                                                                  # index_footprint is sequential

    plt.subplot(2,2,1)
    
#    fig, ax = plt.subplots()

    # Reset the 'sum' image, which is sum of all frames at this footprint.
    # We set this to None rather than a zero array, since we don't know the size of the array until we've 
    # opened the files.
    
    im_sum     = None   

    for index_image in index_images:   # Loop over the images in this observation
        
        t_i = t_group[index_image]  # Grab this, read-only, since we use it a lot.
        
        arr = hbt.read_lorri(t_i['Filename'])
        # arr = hbt.lorri_destripe(arr)

        # Save various parameters in a table, which we will output afterwards
        
        file_arr.append(t_i['Shortname'])
        format_arr.append(t_i['Format'])
        phase_arr.append(t_i['Phase']*hbt.r2d)
        exptime_arr.append(t_i['Exptime'])
        et_arr.append(t_i['ET'])
        utc_arr.append(t_i['UTC'])
        index_footprint_arr.append(index_footprint+1)
        index_hbt_arr.append(f'{index_group}/{index_image}')
        
        if im_sum is None:
            im_sum = arr
        else:    
            im_sum += arr

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = WCS(t_i['Filename'])
    
        # Get the RA / Dec for the four corner points
        
        if t_i['Format'] == '1X1':
            radec_corners = w.wcs_pix2world([[0,0],[0,1024],[1023,1023],[1023,0],[0,0]],0) * hbt.d2r  # Radians
            radec_center  = w.wcs_pix2world([[511,511]],0) * hbt.d2r                       # Radians
        if t_i['Format'] == '4X4':
            radec_corners = w.wcs_pix2world([[0,0],[0,256],[256,256],[256,0],[0,0]],0) * hbt.d2r  # Radians
            radec_center  = w.wcs_pix2world([[128,128]],0) * hbt.d2r                       # Radians
        
    #    ra  = radec[:,0] 
    #    dec = radec[:,1]
    
        # Get RA / Dec for Jupiter
        
        (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'New Horizons')
        (_, ra_jup, dec_jup) = sp.recrad(vec[0:3])  # Radians
        vec_nh_jup = vec[0:3]
         
        plt.plot((radec_corners[:,0]-ra_jup)*hbt.r2d, (radec_corners[:,1]-dec_jup)*hbt.r2d,
                 label = f'{index_group}/{index_image}') # Plot all the four points for this box
                
        is_limb_left = (radec_corners[0,0]-ra_jup)>0
        
        if is_limb_left:
            name_limb_arr.append('Left')
        else:
            name_limb_arr.append('Right')

        # Get the distance from Jupiter center to image center
        
        vec_nh_center = sp.radrec(1, radec_center[0][0], radec_center[0][1])  # Vector from NH, to center of LORRI frame
        ang_jup_center = sp.vsep(vec_nh_jup, vec_nh_center)                   # Ang Sep btwn Jup and LORRI, radians
        dist_jup_center_rj = ang_jup_center * sp.vnorm(vec_nh_jup) / rj_km    # Convert from radians into RJ
        width_lorri = 0.3*hbt.d2r                                             # LORRI full width, radians
    
        dist_jup_center_rj_range = np.array([ang_jup_center+width_lorri/2, ang_jup_center-width_lorri/2]) * \
            sp.vnorm(vec_nh_jup) / rj_km                                      # Finally, we have min and max dist
                                                                              # in the central row of array, in rj
   
        range_rj_arr.append(sp.vnorm(vec_nh_jup)/rj_km)
        
        # If we are on the RHS limb, then we need to reverse this, and negate it.
    
        if not(is_limb_left):
            dist_jup_center_rj_range = -1 * dist_jup_center_rj_range[::-1]

        dist_proj_rj_arr.append((dist_jup_center_rj_range))
        
        imagenum +=1

    dist_proj_rj_foot_arr.append(dist_jup_center_rj_range)  # Set the projected distnace, for the footprint
    name_limb_foot_arr.append(name_limb_arr[-1]) # Defint the limb direction for the footprint
    
    # Finished this loop over image set (aka footprint)
    
    corner_arr = np.array(corner_arr) / rj_km  # This is now an array (n_pts x 5, 3)
    
    # Calculate Jupiter's radius, in degrees
    
    ang_jup_deg = (rj_km / sp.vnorm(vec_nh_jup)) * hbt.r2d
    
    circle1=plt.Circle((0,0),ang_jup_deg,   color='red') # Plot Jupiter itself
    circle2=plt.Circle((0,0),ang_jup_deg*2, facecolor='white', edgecolor='red') # Plot Jupiter itself
    circle3=plt.Circle((0,0),ang_jup_deg*3, facecolor='white', edgecolor='red') # Plot Jupiter itself
    
    plt.gca().add_artist(circle3)
    plt.gca().add_artist(circle2)
    plt.gca().add_artist(circle1)
    plt.gca().set_aspect('equal')
    plt.xlabel('Delta RA from Jupiter  [deg]')    
    plt.ylabel('Delta Dec from Jupiter [deg]')
    plt.xlim((4,-4))
    plt.ylim((-2,2))
    plt.legend(loc='upper right')
    plt.title(f'Gossamer images, {index_group}/{np.amin(index_images)}-{np.amax(index_images)}, f{index_footprint}')
    
    # Set the 'extent', which is the axis values for the X and Y axes for an imshow().
    # We set the y range to be the same as X, to convince python to keep the pixels square.
    
    extent = [dist_jup_center_rj_range[0], dist_jup_center_rj_range[1],\
              dist_jup_center_rj_range[0], dist_jup_center_rj_range[1]]
    
    
#    plt.show() 
    
    # Plot the image, sfit subtracted
    
    plt.subplot(2,2,2)
    degree=5
    im_sum /= len(index_images)  # We have just summed, so now divide
    im_sum_s = hbt.remove_sfit(im_sum,degree=degree)  # _s = sfit subtracted
    plt.imshow(stretch(im_sum_s), origin='lower', extent=extent)
#    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.title(f'{index_group}/{np.amin(index_images)}-{np.amax(index_images)} - p{degree}')

    # Plot the image, with better subtraction
    
    plt.subplot(2,2,3)
    
    method = 'String'
    vars = '6/63-70 p5'
#    index_image = 124
#    index_group = 6
    im_process = hbt.nh_jring_process_image(im_sum, method, vars, 
                                     index_group=index_group, index_image=index_image, mask_sfit=None)


    plt.imshow(stretch(im_process), origin='lower', extent=extent)
#    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.title(f'{index_group}/{np.amin(index_images)}-{np.amax(index_images)} - {vars}')

    im_sum_s_arr.append(im_sum_s)      # Raw summed image from the footprint. Better.
    im_process_arr.append(im_process)  # 'Processed' image from the footprint. Not very good.
        
    # Now load Lauer's processed image of the same footprint
    
    ifp1 = index_footprint+1  # There is an offset of one between my indices, and Tod's
    
    file_lauer = os.path.join(dir_lauer, f'f{ifp1:02}_av_v1f.fits')
    try:
        lun_lauer = fits.open(file_lauer)
        im_lauer = lun_lauer[0].data
        lun_lauer.close()
        plt.subplot(2,2,4)
        plt.imshow(stretch(im_lauer), origin='lower', extent=extent)
        plt.title(f'Lauer f{index_footprint:02}')
        im_lauer_arr.append(im_lauer)
        
    except:   # Skip if, if Lauer did not process it
        im_lauer_arr.append([])

    # Now show all plots

    plt.tight_layout()
    
    plt.show()
    
#    print(f'Is left: {is_limb_left}')
#    print('---')

# Now make the output table for Tod

# =============================================================================
# Make up a table, for Tod Lauer. We don't really used this table -- it's just for output.
# I want to group things     
# =============================================================================

#t_out = Table([name, rho, r, a], names=['Name', 'rho', 'radius', 'a'])
  
t = astropy.table.Table([hbt.frange(1,imagenum), index_hbt_arr, index_footprint_arr, 
                         file_arr, name_limb_arr, dist_proj_rj_arr, phase_arr, 
                         range_rj_arr, et_arr, utc_arr, 
                         format_arr, exptime_arr,],
                        names = ['#', 'HBT #', 'Foot #', 'Name', 'Limb', 'Dist_Proj_RJ', 'Phase', 'Range_RJ',
                                 'ET', 'UTC', 'Format', 'Exptime'])
t['Dist_Proj_RJ'].format='5.3f'
t['Phase'].format = '5.1f'
t['Range_RJ'].format = '4.1f'

path_out = os.path.join(dir_out, 'obstable_gossamer_hbt.txt')
t.write(path_out, format = 'ascii.fixed_width', overwrite='True')
print(f'Wrote: {path_out}')

# =============================================================================
# Now merge all of the footprints into one combined image.
# =============================================================================

# Ideally I would write this as a function. It does not make sense to use WCS for this.
# The coordinate are weird. But it does make sense

# For here, x means 'horizontal axis when plotted', ie, 'second axis of a 2D array'

dx_dist_proj_rj_out = [1, 3]    # Fot consistentency, we should go from -3 to 3.  Negative on left.
dy_dist_proj_rj_out = [0.5, 2]
 
dx_rj        = 0.001          # Resolution of the output image, in RJ. To change output array size, change this.
dy_rj        = dx_rj
dx_pix = int( (dx_dist_proj_rj_out[1] - dx_dist_proj_rj_out[0])/dx_rj ) + 1
dy_pix = int( (dy_dist_proj_rj_out[1] - dy_dist_proj_rj_out[0])/dy_rj ) + 1
vals_x_rj    = hbt.frange(dx_dist_proj_rj_out[0], dx_dist_proj_rj_out[1], dx_pix)
vals_y_rj    = hbt.frange(dy_dist_proj_rj_out[0], dy_dist_proj_rj_out[1], dy_pix)

# Make the stack

stack  = np.zeros((num_footprints, dy_pix, dx_pix))
stack[:,:,:] = np.nan  # Fill it with NaN's

for i in range(num_footprints):
    dist_proj_rj_foot_i = dist_proj_rj_foot_arr[i]  # Get the RJ range for this image (all positive)
    # if ('Right' in name_limb_foot_arr[i]):
    #     dist_proj_rj_foot_i *= -1
    drj_i = np.amax(dist_proj_rj_arr[i]) - np.amin(dist_proj_rj_arr[i])
    fac_zoom = (drj_i / dx_rj) / hbt.sizex_im(im_process_arr[i]) 
    
    # Now do the actual resizing of the image
    
    # im_resize = scipy.ndimage.zoom(im_sum_s_arr[i], fac_zoom)
    
    stretch_percent2 = 95    
    stretch2 = astropy.visualization.PercentileInterval(stretch_percent2) # PI(90) scales to 5th..95th %ile.

    im_resize = scipy.ndimage.zoom(stretch2(im_sum_s_arr[i]), fac_zoom)
 
    DO_OUTPUT_LAUER = False

    if DO_OUTPUT_LAUER:
        
        if len(im_lauer_arr[i]) > 0:
            im_resize = scipy.ndimage.zoom(im_lauer_arr[i], fac_zoom)
        
    # else:
        
    # Now place it on the output stack
    
    x0 = np.digitize(np.amin(dist_proj_rj_foot_i), vals_x_rj) + 0  # add zero to convert from 'array' into int
    x1 = x0 + hbt.sizex_im(im_resize)
    
    y0 = 50 + int(x0/5) # We'll have to figure out a vertical position, but right now we don't.
    
    # y0 = 50 + i*10
    y1 = y0 + hbt.sizey_im(im_resize)

    # print(f'Footprint {i}: pos (y0,x0) = {y0,x0}; shape(large) = {np.shape(large)}; shape(small) = {np.shape(im_resize)}')

    stack[i,:,:] = hbt.replace_subarr_2d(stack[i,:,:], im_resize, (y0,x0))
    # stack[i,:,:] = replace_subarr_2d(stack[i,:,:], im_resize, (y0,x0))
    
    # stack[i, y0:y1, x0:x1] = im_resize
    
# else:
#         print(f'Footprint {i}: No Lauer image')

extent = [ dx_dist_proj_rj_out[0], dx_dist_proj_rj_out[1], dy_dist_proj_rj_out[0], dy_dist_proj_rj_out[1] ]

plt.imshow((np.nanmean(stretch(stack),axis=0)), extent=extent, cmap='plasma', origin='lower')

file_out = os.path.join(dir_out, 'gossamer.png')
plt.savefig(file_out)
print(f'Wrote: {file_out}')
plt.show()


    
#%%%




# And hang on. We've already made backplanes for everything. Can I use those here? 
# I don't know. These images are obviously very edge-on.

#        file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')
#
#        # Save the shortname associated with the current backplane. 
#        # That lets us verify if the backplane for current image is indeed loaded.
#
#        file_backplane_shortname = t_group['Shortname'][index_image]
#				
#        if (os.path.isfile(file_backplane)):
#            lun = open(file_backplane, 'rb')
#            planes = pickle.load(lun)
#            lun.close()

# I could do RA/Dec relative to RA/Dec of Jupiter.
    