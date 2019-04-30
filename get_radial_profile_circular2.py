#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:20:36 2018

@author: throop
"""

import hbt
import numpy as np

def get_radial_profile_circular2(arr, pos = (0,0), width=1, a_xy = (1,1), method='median'):

    """
    Extract a radial profile from an image. 
    
    Will handle a circular case.
    
    With an optional set of semi-major axes passed, will handle a simple elliptical ring.
    
    Parameters
    -----
    
    arr:
        Array of data values (e.g., the image).
        
    pos:
        Position of the center of the circle, in pixels
        
    width:
        Some sort of scaling? Not really sure. Just keep at default.
        
    a_xy:
        Semimajor axes. Tuple. This account for, in a crude way, the tilt of a flat ring.
        e.g., `a_xy = (1, cos(30 deg))` is for a ring tilted so as to be squished vertically.
        Default is for a face-on ring. This will break down for large angles. 
        What matters is the ratio of the supplied x and y values; the actual do not matter.
        
    """

    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    dist_pix_2d = np.sqrt( (((xx - pos[0])/a_xy[1]) ** 2) + (((yy - pos[1])/a_xy[0]) ** 2) )

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


def get_radial_profile_circular_quadrant(arr, pos = (0,0), width=1, a_xy = (1,1), method='median'):

    """
    Extract a radial profile from an image. 
    
    Will handle a circular case.
    
    With an optional set of semi-major axes passed, will handle a simple elliptical ring.
    
    Parameters
    -----
    
    arr:
        Array of data values (e.g., the image).
        
    pos:
        Position of the center of the circle, in pixels
        
    width:
        Some sort of scaling? Not really sure. Just keep at default.
        
    a_xy:
        Semimajor axes. Tuple. This account for, in a crude way, the tilt of a flat ring.
        e.g., `a_xy = (1, cos(30 deg))` is for a ring tilted so as to be squished vertically.
        Default is for a face-on ring. This will break down for large angles. 
        What matters is the ratio of the supplied x and y values; the actual do not matter.
        
    """

    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    dist_pix_2d = np.sqrt( (((xx - pos[0])/a_xy[1]) ** 2) + (((yy - pos[1])/a_xy[0]) ** 2) )

    dist_pix_1d = hbt.frange(0, int(np.amax(dist_pix_2d)/width))*width
        
    profile_1d    = np.zeros((4, len(dist_pix_1d)))
    
    is_q1 = np.logical_and( (xx >= pos[0]), (yy >= pos[1]) )
    is_q2 = np.logical_and( (xx <  pos[0]), (yy >= pos[1]) )
    is_q3 = np.logical_and( (xx <  pos[0]), (yy <  pos[1]) )
    is_q4 = np.logical_and( (xx >= pos[0]), (yy <  pos[1]) )
     
    quadrant_plane = (1 * is_q1) + (2 * is_q2) + (3 * is_q3) + (4 * is_q4)
    
    for i in range(len(dist_pix_1d)-2):
        
        for j in [1,2,3,4]:  # Loop over the quadrants 

    # Identify the pixels which are at the right distance
    
            is_good = np.logical_and( np.logical_and(dist_pix_2d >= dist_pix_1d[i],
                                                   dist_pix_2d <= dist_pix_1d[i+1]),
                                                     quadrant_plane == j)
        
            if (method == 'mean'):
                profile_1d[j-1,i]   = np.nanmean(arr[is_good])
        
            if (method == 'median'):
                profile_1d[j-1,i] = np.nanmedian(arr[is_good])
            
    return (dist_pix_1d, profile_1d)

if (__name__ == '__main__'):
    
    # Make a circular profile
    
    arr = hbt.dist_center(100)
    
    (radius, profile) = get_radial_profile_circular(arr, pos=(50,50))
    plt.plot(radius,profile)
    plt.xlabel('Radius')
    plt.ylabel('Value')
    plt.show()
    
    # Make a quadrant profile
    
    (radius, profile) = get_radial_profile_circular_quadrant(arr, pos=(45,50))
    for i in range(4):
        plt.plot(radius,profile[i,:])
        
    plt.xtitle('Radius')
    plt.ytitle('Value')
    plt.show()
    
    
    
    pass

