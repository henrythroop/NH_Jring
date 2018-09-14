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

# =============================================================================
# This is the main Gossamer ring image processing / mosaic code. 
# It reads and processes raw FITS files and outputs a merged FITS file mosaic.
# It also reads and processes the PCA-processed stacks of Gossamer images that Tod Lauer has made
#
# September 2018 HBT
# =============================================================================
# 
def find_file_footprint_lauer(dir_lauer, str_met):
    
    """
    Given a string for an MET, finds which of Tod's footprints corresponds to it.
    This is a simple one-off function. A single line of 'grep' would do the same, but I couldn't get it to work.
    
    Tod made a bunch of PCA'd images for me, which are in FITS files like 'f02_av_v3.fits'. The contents of tehse
    files are in accompanying 'f03.lis' text file. This function here searches for matches.
    
    """
    
    files = glob.glob(os.path.join(dir_lauer, 'f[0-2]*.lis'))  # Tod made a bunch of text files, listing the MET's.
    for file in files:
        
        if str_met in open(file).read():  # This is great! So pythonic and easy.
            file_fits = file.replace('.lis', '_av_v3.fits')   # Create the filename of the output stack.
            return file_fits
    
    return None

#%%%
    
# NB: This program is just code. It is not one callable function. That's OK, and makes debugging easier, I guess.
    
plt.set_cmap('plasma')

rj_km = 71492

# List all of the Gossamer observations. I have determined these groupings manually, based on timings.
# In general it looks like usually these are stacks of four frames at each pointing ('footprint')

# Lauer's footprints are just numbered sequentially (f01, f02, etc). So, we need to keep these files here in 
# same order, so as to correlate exactly w/ Lauer's.

index_group = 6

index_image_list = [
                    np.array([54]),    # 1X1. Main ring. 1.6 .. 1.9 RJ. Right ansa. One image, edge-on.
                    hbt.frange(59,62),  # Right ansa. -2.7 RJ. Gossamer maybe there but hidden -- need to optimize.
                    hbt.frange(63,66),  # Left ansa. Furthest out, no ring. Use for subtraction.
                    hbt.frange(67,70),  # Left ansa. Far out, no ring. Closer than above.
                    hbt.frange(71,74),  # Left ansa. No ring. Closer than above.
                    hbt.frange(75,78),  # Left ansa. Gossamer. Closer than above.
                    np.array([80]),  # Left ansa. Main ring. Prob gossamer too but not interested. 1X1. Svrl ptgs.
                    hbt.frange(95,98),  # Right ansa. Gossamer limb.
                    hbt.frange(99,102), # Left ansa. Similar geometry as 75-78.
                    hbt.frange(112,115),# Right ansa. Gossamer limb. Similar geometry as 95-98 (offset pointing)
                    hbt.frange(116,119),# Left ansa. Similar geometry as 75-78. Sat visible??
                    hbt.frange(120,123),  # Left ansa. Closer than above.
                    hbt.frange(124,127),  # Left ansa. Closer than above.
                    hbt.frange(128,131),  # Left ansa + main ring ansa. Closer than above.
                    hbt.frange(132,135),  # Main ring. Left side, well inside of gossamer.
                    
                    hbt.frange(157,160),
                    hbt.frange(173,176),
                    hbt.frange(177,180),
                    hbt.frange(181,184),
                    hbt.frange(185,188),
                    hbt.frange(207,210),
                    hbt.frange(211,214),
                    hbt.frange(225,227),
                    np.array([228]),
                    hbt.frange(229,232),   # Right gossamer ansa.
                    # hbt.frange(233,238),   # Main ring? Many images, lots of jitter -- can't stack in place.
                    
                    ]
                    

# index_image_list = [
#                       # hbt.frange(46,48), # 1X1. Main ring. Not edge-on, not part of this sequence.
#                       # hbt.frange(49,53), # 1X1. Main ring. Not edge-on, not part of this sequence.
#                     hbt.frange(79,88),  # Left ansa. Main ring. Not edge-on. 1X1. Svrl ptgs.
                       # np.array([54]),    # 1X1. Main ring. 1.6 .. 1.9 RJ. Right ansa. One image, edge-on.
                       # np.array([89]),    # 1X1. Main ring. 1.6 .. 1.9 RJ. Right ansa. One image, edge-on.
#                       # hbt.frange(59,62),  # Right ansa. 2.4 .. 2.7 RJ. Gossamer there, bit hidden?
# #                       hbt.frange(63,66),  # Left ansa. -3.2 .. 3.0. Furthest out, no ring. Use for subtraction.
# #                       hbt.frange(67,70),  # Left ansa. -3.0 .. -2.8. Far out, no ring.
# #                       hbt.frange(71,74),  # Left ansa. -2.8 .. -2.6. Far out, no ring.
# #                       hbt.frange(75,78),  # Left ansa. Gossamer. -2.6 .. -2.4 RJ. Closer than above.
#                     hbt.frange(89,94),  # Io! 1.6-1.85. Right ansa. 1X1. Lots jitter. Could do one frame only?

# #                       hbt.frange(59,62),  # Right ansa.           2.4 .. 2.7. Gossamer is there but hidden?
# #                       hbt.frange(95,98),  # Right ansa. Gossamer limb. 2.35 .. 2.65. Definitely gossamer.
# #                       hbt.frange(99,102), # Left ansa. Gossamer. 2.65 .. 2.4. 
# #                       hbt.frange(112,115),# Right ansa. Gossamer limb. 2.4 .. 2.6. Def Gossamer.
# #                        hbt.frange(116,119),  # Left ansa. -2.65 .. -2.35
# #                        hbt.frange(120,123),  # Left ansa. -2.45 .. -2.15. 
# #                        hbt.frange(124,127),  # Left ansa. -2.25 .. -1.95. 
# #                        hbt.frange(128,131),  # Left ansa. -2 .. -1.70. Gossamer + main.
                    
# # # #                      hbt.frange(157,160),# Left ansa.           -2.7 .. 2.4
# # # #                     hbt.frange(173,176),# Right ansa. Main ring. 1.8 .. 2.
# # # #                      hbt.frange(177,180),# Right ansa. Gossamer.  2.0 .. 2.3.
# # # #                      hbt.frange(181,184), # Right ansa. Gossamer. 2.2 .. 2.5.
# # # #                      hbt.frange(185,188),   # Right ansa. Gossamer. 2.4 .. 2.8. Outer edge of Gossamer?
# # # #                      hbt.frange(207,210), # Right ansa. Gossamer. 2.35 .. 2.65
# # # #                      hbt.frange(211,214),  # Left ansa. -2.65 .. -2.40. Outer edge of Gossamer?
# # # #                      hbt.frange(225,227),    # Left ansa. -2.7 .. 2.35. Outer edge of Gossamer?
# # # #                      np.array([228]),           #  Left ansa. -2.7 .. 2.35. Outer edge of Gossamer?
# # # # #                     # hbt.frange(229,232),     # Right ansa. 2.45-2.70. Outer edge of Gossamer.
# # # # #                     # hbt.frange(233,238),   # 1.6 .. 1.9. Main ring. Probably should not plot these -- dfft sequence.
# # # # # #                    
                    # ]
                    


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

ra_arr           = []
dec_arr          = []
dist_proj_rj_arr = []      # Projected distance in RJ, for each image
dist_proj_rj_foot_arr = [] # Projected distance in RJ, for each footprint (ie, set of four images)
dist_jup_center_vertical_rj_foot_arr = []
dist_jup_center_horizontal_rj_foot_arr = []
range_rj_arr     = []
angle_elev_arr   = []
et_arr           = []
dx_arr           = []
dy_arr           = []
dz_arr           = []
corner_arr       = []
format_arr       = []
exptime_arr      = []
name_limb_arr    = []
name_limb_foot_arr=[]
file_arr         = []
phase_arr        = []
utc_arr          = []
index_hbt_arr    = []
index_footprint_arr=[]
im_process_arr   = []   # Array of final processed images
im_lauer_arr     = []
im_sum_s_arr     = []
    
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
            
            arr *= 16
            
        if t_i['Format'] == '4X4':
            radec_corners = w.wcs_pix2world([[0,0],[0,256],[256,256],[256,0],[0,0]],0) * hbt.d2r  # Radians
            radec_center  = w.wcs_pix2world([[128,128]],0) * hbt.d2r                       # Radians
    
        # Get RA / Dec for Jupiter
        
        (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'New Horizons')
        (_, ra_jup, dec_jup) = sp.recrad(vec[0:3])  # Radians
        vec_nh_jup = vec[0:3]
        
        # Get vectors Sun-Jup, and Sun-NH
        
        (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'Sun')
        vec_sun_jup = vec[0:3]  # Radians
        
        vec_sun_nh = vec_sun_jup - vec_nh_jup
        
        # Plot the boxes
        
        plt.plot((radec_corners[:,0]-ra_jup)*hbt.r2d * (-1), (radec_corners[:,1]-dec_jup)*hbt.r2d,
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
    
        dist_jup_center_rj_range = np.array([ang_jup_center-width_lorri/2, ang_jup_center+width_lorri/2]) * \
            sp.vnorm(vec_nh_jup) / rj_km                                      # Finally, we have min and max dist
                                                                              # in the central row of array, in rj
        range_rj_arr.append(sp.vnorm(vec_nh_jup)/rj_km)
   
        # Calc the elevation angle (aka sub-obs lat on Jup)
        
        mx = sp.pxform('J2000', 'IAU_JUPITER', t_i['ET'])
        vec_nh_jup_jup = sp.mxv(mx, vec_nh_jup)
        (dist, lon, lat) = sp.recrad(-vec_nh_jup_jup)
        
        angle_elev_arr.append(lat * hbt.r2d)  # Save the sub-obs latitude, in degrees
        
        # If we are on the LHS limb, then we need to reverse this, and negate it.
    
        if (is_limb_left):
            dist_jup_center_rj_range = -1 * dist_jup_center_rj_range[::-1]

        dist_proj_rj_arr.append((dist_jup_center_rj_range))
        
        # Get the distance above the ring plane -- that is, distance from image center to ring plane, projected.
        # To do this:
        #   - Project a vector from s/c along central LORRI pixel.
        #   - Make a SPICE 'line' from s/c along this vector
        #   - Get C/A point from this 'line', to Jupiter center
        #   - Convert that point to a distance and elevation using RECRAD, in IAU_JUPITER coords.
        
        # Get the 'close point', which is basically the point in space above/below the ring plane, where NH is pointed.
        # "Finds the nearest point on a line." So I guess it is in J2000 from the sun position, baiscally.
        
        (pt_closest, dist) = sp.nplnpt(vec_sun_nh, vec_nh_center, vec_sun_jup)
        
        # Convert this to a vector centered on Jupiter
        
        pt_closest_jup = -pt_closest + vec_sun_jup
        
        # Now convert this to an elevation angle and distance
        # 'Radius' will be the centerpoint of the image, in km from Jup.
        
        mx = sp.pxform('J2000', 'IAU_JUPITER', t_i['ET'])
        
        # Convert to IAU_JUP coords. Looks good.
        
        pt_closest_jup_jup = sp.mxv(mx, pt_closest_jup)  # Get the close point to Jup, in IAU_JUP coords
        
        # Convert into radius / lon / lat, in km, in Jupiter frame. 
        # ** Seems to be some error here. Radius is OK, but lat not right??
        
        (radius, lon, lat) = sp.reclat(pt_closest_jup_jup)
        
        dist_jup_center_horizontal_rj = math.cos(lat) * radius / rj_km   # This is now the same as dst_jup_center_rj_range. Looks right.
        dist_jup_center_vertical_rj   = math.sin(lat) * radius / rj_km
        
        imagenum +=1 # Loop to next image number in the footprint
        
    dist_proj_rj_foot_arr.append(dist_jup_center_rj_range)  # Set the projected distance, for the footprint
    name_limb_foot_arr.append(name_limb_arr[-1])            # Define the limb direction for the footprint
    
    dist_jup_center_horizontal_rj_foot_arr.append(dist_jup_center_horizontal_rj)
    dist_jup_center_vertical_rj_foot_arr.append(dist_jup_center_vertical_rj)
    
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
    plt.xlim((-4,4))
    plt.ylim((-2,2))
    plt.legend(loc='upper right')
    plt.title(f'Gossamer images, {index_group}/{np.amin(index_images)}-{np.amax(index_images)}, f{index_footprint}')

    print(f'Position of LORRI center in Jup coords: {pt_closest_jup_jup}')
    print(f'Position of LORRI center: Horiz = {dist_jup_center_horizontal_rj:.2}, Vert = {dist_jup_center_vertical_rj:.2}')
   
    # Set the 'extent', which is the axis values for the X and Y axes for an imshow().
    # We set the y range to be the same as X, to convince python to keep the pixels square.
    
    extent = [dist_jup_center_rj_range[0], dist_jup_center_rj_range[1],\
              dist_jup_center_rj_range[0], dist_jup_center_rj_range[1]]
        
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
    
    files_short = t_group[index_image_list[index_footprint]]['Shortname']
    file_short = files_short[0]
    met = file_short.split('_')[1]  # Get the MET embedded in the 0th filename in this footprint
    
    file_lauer = find_file_footprint_lauer(dir_lauer, met)
    
    if file_lauer:
        lun_lauer = fits.open(file_lauer)
        im_lauer = lun_lauer[0].data    # Load Lauer's footprint
        lun_lauer.close()
        plt.subplot(2,2,4)
        plt.imshow(stretch(im_lauer), origin='lower', extent=extent)
        plt.title(f'Lauer f{index_footprint:02}')
        im_lauer_arr.append(im_lauer)
        
    else:   # Skip if, if Lauer did not process it
        im_lauer_arr.append([])

    # Now show all plots

    plt.tight_layout()
    
    plt.show()
    
#    print(f'Is left: {is_limb_left}')
#    print('---')


#%%%

# =============================================================================
# Now merge all of the footprints into one combined image.
# =============================================================================

# Ideally I would write this as a function. It does not make sense to use WCS for this.
# The coordinate are weird. But it does make sense

# For here, x means 'horizontal axis when plotted', ie, 'second axis of a 2D array'

do_all = False
do_rhs_only = False
do_lhs_only = True

if (do_lhs_only):
    dx_dist_proj_rj_out = [-3, 3]    # Hortizontal size of output plot, in RJ
    dy_dist_proj_rj_out = [-0.5, 0.5]    # Vertical    size of output plot, in RJ
 
dx_rj        = 0.001          # Resolution of the output image, in RJ. To change output array size, change this.
dy_rj        = dx_rj
dx_pix = int( (dx_dist_proj_rj_out[1] - dx_dist_proj_rj_out[0])/dx_rj ) + 1
dy_pix = int( (dy_dist_proj_rj_out[1] - dy_dist_proj_rj_out[0])/dy_rj ) + 1
vals_x_rj    = hbt.frange(dx_dist_proj_rj_out[0], dx_dist_proj_rj_out[1], dx_pix)
vals_y_rj    = hbt.frange(dy_dist_proj_rj_out[0], dy_dist_proj_rj_out[1], dy_pix)

# =============================================================================
# Make the stack, combining the Throop and Lauer mosaics into one array
# =============================================================================

stack  = np.zeros((num_footprints, dy_pix, dx_pix))
stack[:,:,:] = np.nan  # Fill it with NaN's

stack_lauer = stack.copy()

for i in range(num_footprints):

    dist_proj_rj_foot_i = dist_proj_rj_foot_arr[i]  # Get the RJ range for this image (all positive)

    drj_i = np.amax(dist_proj_rj_arr[i]) - np.amin(dist_proj_rj_arr[i])
    fac_zoom = (drj_i / dx_rj) / hbt.sizex_im(im_process_arr[i]) 
    
    # Resize the image. Ignore a deprecation warning.
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")    
        im_resize = scipy.ndimage.zoom(im_process_arr[i], fac_zoom)

# Now place it on the output stack

    x0 = np.digitize(np.amin(dist_proj_rj_foot_i), vals_x_rj) + 0  # add zero to convert from 'array' into int
    x1 = x0 + hbt.sizex_im(im_resize)
    
    y0 = int(np.digitize(-dist_jup_center_vertical_rj_foot_arr[i], vals_y_rj) - hbt.sizey(im_resize)/2  + 0)
    y1 = y0 + hbt.sizey_im(im_resize)
    
    stack[i,:,:] = hbt.replace_subarr_2d(stack[i,:,:], im_resize, (y0,x0))

    # print(f'Now placing footprint {i} at location x0={x0}, y0={y0} in box of size {np.shape(stack)}')

# Now place the Lauer image, if it exists
        
    if len(im_lauer_arr[i]) > 0:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            im_resize_lauer = scipy.ndimage.zoom(im_lauer_arr[i], fac_zoom)
                    
        # Place it on the output stack
        
        stack_lauer[i,:,:] = hbt.replace_subarr_2d(stack_lauer[i,:,:], im_resize_lauer, (y0,x0))
                     
    extent = [ dx_dist_proj_rj_out[0], dx_dist_proj_rj_out[1], dy_dist_proj_rj_out[0], dy_dist_proj_rj_out[1] ]

# Flatten the stack, using mean() or median()
    
im       = np.nanmean(stack,       axis=0)
im_lauer = np.nanmean(stack_lauer, axis=0)

# Plot the flattened stacks, if requested

do_plot_stacks = False

if do_plot_stacks:
    plt.imshow( stretch(im), extent=extent, cmap='plasma', origin='lower')
    plt.title('Throop Processed')
    file_out = os.path.join(dir_out, 'gossamer_processed.png')
    plt.savefig(file_out)
    print(f'Wrote: {file_out}')
    plt.show()
        
    do_plot_lauer = True
    
    if do_plot_lauer:
        if np.nansum(stack_lauer) > 0:  # If there is some data in it
            plt.imshow(stretch(im_lauer), extent=extent, cmap='plasma', origin='lower')
            plt.title('Lauer')
            file_out = os.path.join(dir_out, 'gossamer_lauer.png')
            plt.savefig(file_out)
            print(f'Wrote: {file_out}')
            plt.show()
        
# Now take these images and re-assemble them into one mosaic of Lauer + Throop together
            
im_mosaic = np.vstack((im, im_lauer))

# Now collapse the whitespace in the merged image, by removing any empty rows and column

im_mosaic[0:2,2000:2100] = 0.01 # Fill in a little region to mask where it will not be collapsed

# Find the empty rows and columns, and remove them.

rows = np.nansum(im_mosaic, axis=1)
columns = np.nansum(im_mosaic, axis=0)
im_mosaic_crop = np.delete(im_mosaic,      np.where(rows==0)[0], axis=0)
im_mosaic_crop = np.delete(im_mosaic_crop, np.where(columns==0)[0], axis=1)

# Plot the image -- normal

hbt.figsize((15,15))

plt.imshow(im_mosaic_crop, origin='lower')
plt.show()

# Plot the image -- stretched

plt.imshow(stretch(im_mosaic_crop), origin='lower')
plt.show()

hbt.figsize()

# And save to disk as a FITS file

file_out = 'gossamer_mosaic.fits'
hdu = fits.PrimaryHDU(im_mosaic_crop)
hdu.writeto(os.path.join(dir_out, file_out), overwrite=True)
print(f'Wrote: {os.path.join(dir_out,file_out)}')

#%%%

# =============================================================================
# Make up a table, for Tod Lauer. We don't really used this table -- it's just for output.
# I want to group things     
# =============================================================================

#t_out = Table([name, rho, r, a], names=['Name', 'rho', 'radius', 'a'])
  
t = astropy.table.Table([hbt.frange(1,imagenum), index_hbt_arr, index_footprint_arr, 
                         file_arr, name_limb_arr, dist_proj_rj_arr, phase_arr, 
                         range_rj_arr, angle_elev_arr, et_arr, utc_arr, 
                         format_arr, exptime_arr,],
                        names = ['#', 'HBT #', 'Foot #', 'Name', 'Limb', 'Dist_Proj_RJ', 'Phase', 'Range_RJ',
                                 'Elev', 'ET', 'UTC', 'Format', 'Exptime'])
t['Dist_Proj_RJ'].format='5.3f'
t['Phase'].format = '5.1f'
t['Elev'].format = '5.1f'
t['Range_RJ'].format = '4.1f'

path_out = os.path.join(dir_out, 'obstable_gossamer_hbt.txt')
t.write(path_out, format = 'ascii.fixed_width', overwrite='True')
print(f'Wrote: {path_out}')

#%%%

    