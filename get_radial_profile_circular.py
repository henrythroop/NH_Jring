#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:20:36 2018

@author: throop
"""

import hbt
import numpy as np
import math
import matplotlib.pyplot as plt # pyplot


import glob
import math
import os.path
import os
from astropy.table import Table
import astropy
from   astropy.io import fits
import astropy.table
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library
import pickle # For load/save

from   astropy.wcs import WCS

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.


def get_profile_linear(arr, pos=None, method='median', axis_move = 1):

    """
    Get the linear profile of an array. That is, the rowsum or the column-sum.
    
    Also:
        - Returns the row numbers (or column numbers).
        - Allows summing along either axis.
        - Properly handles case of odd # vs. even # of rows.
    

    Arguments
    -----
    
    arr:
        An array. Typically 2d NumPy array.
    
    Optional keyword arguments:
    -----
    
    pos:
        Positon of the center. Otherwise, it is calculated to be at the middle.
        
    method:
        String. Either 'mean' or 'median'.
        
    axis_move:
        Which axis do we use (ie, move along). 0 or 1, typically. **Needs to be validated -- might be swapped**
        
    """
        
    if not(pos):
        pos = np.shape(arr)[axis_move]/2
        
    dist_pix_1d = np.array(hbt.frange(0, np.shape(arr)[axis_move]-1)).astype(float)
    dist_pix_1d = dist_pix_1d - pos + 0.5

    if method=='median':
        profile = np.nanmedian(arr, axis=axis_move)
    if method=='mean':
        profile = np.nanmean(arr, axis=axis_move)

    # plt.plot(dist_pix_1d, profile)
    # plt.show()
    
    return(dist_pix_1d, profile)    

# =============================================================================
# GET RADIAL PROFILE CIRCULAR 
# =============================================================================
    
def get_radial_profile_circular(arr, pos = None, width=1, a_xy = (1,1), theta = 0, method='median'):

    """
    Extract a radial profile from an image. 
    
    Will handle a circular case.
    
    With an optional set of semi-major axes passed, will handle a simple elliptical ring.
    
    Parameters
    -----
    
    arr:
        Array of data values (e.g., the image).
    
    Optional Keyword Paramters:
        
    pos:
        Position of the center of the circle, in pixels. If omitted, assume center of image.
        
    width:
        Bin width, in pixels. If this is the default, then the length of the output array
        will be essentially the distance from center to diagonal corner -- i.e,., sqrt(2) * hbt.sizex(arr)/2
        
    a_xy:
        Semimajor axes. Tuple. This account for, in a crude way, the tilt of a flat ring.
        e.g., `a_xy = (1, cos(30 deg))` is for a ring tilted so as to be squished vertically.
        Default is for a face-on ring. This will break down for large angles. 
        What matters is the ratio of the supplied x and y values; the actual do not matter.
    
    theta:
        Rotational angle, in radians.
        
    """

    if not pos:
        pos = np.array(np.shape(arr))/2
        
    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    # dist_pix_2d = np.sqrt( (((xx - pos[0])/a_xy[1]) ** 2) + (((yy - pos[1])/a_xy[0]) ** 2) )
    
    dist_pix_2d = np.sqrt(    ( (( ((xx-pos[0])*math.cos(theta) - (yy-pos[1])*math.sin(theta))) / a_xy[1]) ** 2) + 
                              ( (( ((yy-pos[1])*math.cos(theta) + (xx-pos[0])*math.sin(theta))) / a_xy[0]) ** 2)   )
  

    dist_pix_1d = hbt.frange(0, int(np.amax(dist_pix_2d)/width))*width
    
    profile_1d    = 0. * dist_pix_1d.copy()
        
    for i in range(len(dist_pix_1d)-2):

    # Identify the pixels which are at the right distance
    
        is_good = np.logical_and(dist_pix_2d >= dist_pix_1d[i],
                                 dist_pix_2d <= dist_pix_1d[i+1]) 
    
        if (method == 'mean'):
            profile_1d[i]   = np.nanmean(arr[is_good])
    
        if (method == 'median'):
            profile_1d[i] = np.nanmedian(arr[is_good])
        
    return (dist_pix_1d, profile_1d)

# =============================================================================
# GET RADIAL PROFILE CIRCULAR 
# =============================================================================
    
def get_radial_profile_circular_quadrant(arr, pos = None, width=1, a_xy = (1,1), theta = 0, method='median'):

    """
    Extract a radial profile from an image. 
    
    Will handle a circular case.
    
    With an optional set of semi-major axes passed, will handle a simple elliptical ring.
    
    Parameters
    -----
    
    arr:
        Array of data values (e.g., the image).
    
    Optional Keyword Paramters:
        
    pos:
        Position of the center of the circle, in pixels. If omitted, assume center of image.
        
    width:
        Bin width, in pixels. If this is the default, then the length of the output array
        will be essentially the distance from center to diagonal corner -- i.e,., sqrt(2) * hbt.sizex(arr)/2
        
    a_xy:
        Semimajor axes. Tuple. This account for, in a crude way, the tilt of a flat ring.
        e.g., `a_xy = (1, cos(30 deg))` is for a ring tilted so as to be squished vertically.
        Default is for a face-on ring. This will break down for large angles. 
        What matters is the ratio of the supplied x and y values; the actual do not matter.
    
    theta:
        Rotation angle.
        
    """

    if not pos:
        pos = np.array(np.shape(arr))/2
        
    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    # Make the 2D distance array
    
    # Equation for an ellipse, non-rotated
    
    # dist_pix_2d = np.sqrt( (((xx - pos[0])/a_xy[1]) ** 2) + (((yy - pos[1])/a_xy[0]) ** 2) )
    
    # Equation for an ellipse, rotated. This seems to work.
    
    dist_pix_2d = np.sqrt(    ( (( ((xx-pos[0])*math.cos(theta) - (yy-pos[1])*math.sin(theta))) / a_xy[1]) ** 2) + 
                              ( (( ((yy-pos[1])*math.cos(theta) + (xx-pos[0])*math.sin(theta))) / a_xy[0]) ** 2)   )

    # Plot the 2D distance array. This is just for debugging really.
    
    do_plot_2d = False

    if do_plot_2d:
        plt.contour(dist_pix_2d)
        plt.gca().set_aspect('equal')
        plt.show()

    # Make a 1D array, which is a list of distances from 0 .. (max of 2d array)
    
    dist_pix_1d = hbt.frange(0, int(np.amax(dist_pix_2d)/width))*width
        
    profile_1d    = np.zeros((4, len(dist_pix_1d)))
    
    # Define some quadrants. These can really be defined however we like. In this case it is based on
    # non-rotated positions -- just a straightforward xy.
    # (Definitions seems backwards, but they give the right result like this.)
    
    is_q1 = np.logical_and( (xx >= pos[0]), (yy >= pos[1]) )
    is_q4 = np.logical_and( (xx <  pos[0]), (yy >= pos[1]) )
    is_q3 = np.logical_and( (xx <  pos[0]), (yy <  pos[1]) )
    is_q2 = np.logical_and( (xx >= pos[0]), (yy <  pos[1]) )
    
    quadrant_plane = (1 * is_q1) + (2 * is_q2) + (3 * is_q3) + (4 * is_q4)

    # Do an implant, if requested.
    # This is super-simple -- just set to zero the pixels between a set of radii.  

    do_implant = False
    
    if do_implant:
        pixscale_km  = 1.765  # This is pix scale for Teacup image of Tod
        a_implant_km = np.array([1000, 1050])
        a_implant_pix = a_implant_km / pixscale_km

        mask = np.logical_not(np.logical_and(dist_pix_2d > a_implant_pix[0], 
                                             dist_pix_2d < a_implant_pix[1]))
    
        arr = mask * arr
            
    for i in range(len(dist_pix_1d)-2):
        for j in [1,2,3,4]:

    # Identify the pixels which are at the right distance.
    # Loop over the 1D array, and for each distance, find pixels in the 2D array at that distance. 
    
            is_good = np.logical_and( np.logical_and(dist_pix_2d >= dist_pix_1d[i],
                                                     dist_pix_2d <= dist_pix_1d[i+1]), quadrant_plane == j) 
        
            if (method == 'mean'):
                profile_1d[j-1, i]   = np.nanmean(arr[is_good])
        
            if (method == 'median'):
                profile_1d[j-1, i] = np.nanmedian(arr[is_good])
        
    # If requested, make a plot of the quadrants, superimposed with the data array
    
    do_plot_quadrants = False
    
    if do_plot_quadrants:
        plt.imshow(quadrant_plane, origin='lower')
        plt.show()
        plt.imshow(stretch(quadrant_plane * arr), origin='lower')
        plt.show()
        plt.imshow(stretch(arr), origin='lower')
        plt.show()

    # If requested, make a plot of the distance, masked into a ring -- essentially a ring implant
    
    do_plot_quadrants_plus_data = False
    
    if do_plot_quadrants_plus_data and do_implant:
        # plt.imshow(mask, origin='lower')
        # plt.show()
        plt.imshow(stretch(mask * arr), origin='lower')
        plt.show()


    # Return results
    
    return (dist_pix_1d, profile_1d)

# =============================================================================
# End of Function
# =============================================================================

if (__name__ == '__main__'):
    
    # Make a circular profile
    
    arr = hbt.dist_center(101)
    
    (radius, profile) = get_radial_profile_circular(arr, pos=(50,50))
    plt.plot(radius,profile)
    plt.xlabel('Radius')
    plt.ylabel('Value')
    plt.show()
    
    # Make a quadrant profile
    
    (radius, profile) = get_radial_profile_circular_quadrant(arr, pos=(59,52))
    for i in range(4):
        plt.plot(radius,profile[i,:], linewidth=5, alpha=0.3)
        
    plt.xlabel('Radius')
    plt.ylabel('Value')
    plt.show()
        
    pass
