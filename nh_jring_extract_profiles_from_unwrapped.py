#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:41:16 2017

Extracts radial and azimuthal profiles from unwrapped ring images.

bins_radius:  Input: the definition of radius  bins in the unwrapped image
bins_azimuth: Input: The definition of azimuth bins in the unwrapped image

radius_out:   Fraction of radius  range to use, when extracting azimuthal profile [sic]
azimuth_out:  Fraction of azimuth range to use, when extracting radial    profile [sic]

Or, if these are single floating values (e.g., 0.5), then use the central fraction of the pixels available
(e.g., inner 50%)


@author: throop
"""

def nh_jring_extract_profiles_from_unwrapped(im_unwrapped, bins_radius, bins_azimuth, radius_out, azimuth_out, 
                                             mask_unwrapped=False):

    import hbt
    import pickle
    from   astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np
    import os.path
    import astropy
    from   scipy.interpolate import griddata
    import math
    import spiceypy as sp
    import warnings

    if type(mask_unwrapped) == type(np.array([])):
        DO_MASK = True

    print("NH_J_extract_profiles_from_unwrapped: DO_MASK = {}".format(DO_MASK))
     
    if (DO_MASK):    
        is_good_unwrapped = (mask_unwrapped == False)

        im_unwrapped_masked = im_unwrapped.copy()
        im_unwrapped_masked[is_good_unwrapped == False] = np.nan
    
        im2 = im_unwrapped_masked

    else:
        im2 = im_unwrapped
        
# First extract the radial profile. For this, we use a subset of the azimuth angles.
    
    bin_0 = int(hbt.sizex(bins_azimuth) * (0.5 - (0.5 * azimuth_out)))
    bin_1 = int(hbt.sizex(bins_azimuth) * (0.5 + (0.5 * azimuth_out)))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)    
        profile_radius  = np.nanmean(im2[:, bin_0:bin_1],1)
        
    bins_radius_out = bins_radius[bin_0:bin_1]
    
# Now extract the azimuthal profile
    
    bin_0 = int(hbt.sizex(bins_radius) * (0.5 - (0.5 * radius_out)))
    bin_1 = int(hbt.sizex(bins_radius) * (0.5 + (0.5 * radius_out)))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)    
        profile_azimuth = np.nanmean(im2[bin_0:bin_1, :],0)
        
    bins_azimuth_out  = bins_azimuth[bin_0:bin_1]

    plt.plot(bins_azimuth, profile_azimuth)
    plt.title("p_from_unwrapped, Az")
    plt.show()

    plt.plot(bins_radius, profile_radius)
    plt.title("p_from_unwrapped, Rad")
    plt.show()

    
# Now return everything

    return (profile_radius, profile_azimuth)
    