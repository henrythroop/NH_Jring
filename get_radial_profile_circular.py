#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:20:36 2018

@author: throop
"""

import hbt
import numpy as np

def get_radial_profile_circular(arr, pos = (0,0), width=1, a_xy = (1,1), method='median'):

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
            profile_1d[i]   = np.mean(arr[is_good])
    
        if (method == 'median'):
            profile_1d[i] = np.median(arr[is_good])
        
    return (dist_pix_1d, profile_1d)