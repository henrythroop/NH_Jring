#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 11:53:55 2018

@author: throop
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:20:36 2018

@author: throop
"""

import hbt
import numpy as np
import astropy.stats
import warnings

def get_radial_profile_backplane(im, radius_plane, method='median', num_pts = 100, do_std=False):

    """
    Extract a radial profile from an image. 
    
    Uses a backplane passed in.
    
    Parameters
    -----
    
    im:
        Array of data values (e.g., the image).
        
    radius_plane:
        2D array, which is the backplane. Typically this is planes['Radius_eq'].
    
    num_pts:
        Scalar. Number of points to use in the output array.  Output radius is evenly spaced
        from 0 .. max(radius_plane).
        
    Optional parameters
    -----    
    
    method: 
        String. 'mean' or 'median'.
        
    do_stdev: 
        Boolean. If set, compute the standard deviation, and return in the tuple
        
    """
    radius_1d = hbt.frange(0, int(np.amax(radius_plane)), num_pts)
    
    profile_1d    = 0. * radius_1d.copy()
    
    std_1d        = profile_1d.copy()
        
    for i in range(len(profile_1d)-2):

    # Identify the pixels which are at the right distance
    
        is_good = np.logical_and(radius_plane >= radius_1d[i],
                                 radius_plane <= radius_1d[i+1]) 
    
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)    

            if (method == 'mean'):
                profile_1d[i]   = hbt.nanmean(im[is_good])
        
            if (method == 'median'):
                profile_1d[i] = hbt.nanmedian(im[is_good])
    
            # clipped = astropy.stats.sigma_clip(im[is_good], sigma=1.5)
            
            # std_1d[i] = np.nanstd(clipped)
            
    if do_std: 
        return (radius_1d, profile_1d, std_1d)
    
    return (radius_1d, profile_1d)


def get_radial_profile_backplane_quadrant(im, radius_plane, longitude_plane, method='median', num_pts = 100, do_std=False):

    """
    Extract a radial profile from an image. 
    
    Uses a backplane passed in.
    
    Returns four values, corresponding to the four quadrants.
    
    Parameters
    -----
    
    im:
        Array of data values (e.g., the image).
        
    radius_plane:
        2D array, which is the backplane. Typically this is planes['Radius_eq'].
    
    num_pts:
        Scalar. Number of points to use in the output array.  Output radius is evenly spaced
        from 0 .. max(radius_plane).
        
    Optional parameters
    -----    
    
    method: 
        String. 'mean' or 'median'.
        
    do_stdev: 
        Boolean. If set, compute the standard deviation, and return in the tuple
        
    """
    
    radius_1d = hbt.frange(0, int(np.amax(radius_plane)), num_pts)
    
    profile_1d    = np.zeros((4, len(radius_1d)))
    
    std_1d        = profile_1d.copy()
    
    quadrant_plane = np.ceil((np.pi + longitude_plane) * 4 / (np.pi * 2)).astype(int) - 1   # Range = 0 .. 3
    
    for i in range(len(radius_1d)-2):
        
        for j in range(4):    

        # Identify the pixels which are at the right distance
        
            is_good = np.logical_and(np.logical_and(radius_plane >= radius_1d[i],
                                     radius_plane <= radius_1d[i+1],),
                                     quadrant_plane == j)
        
            # print(f'Summing {np.sum(is_good)} cells')
            
            if (method == 'mean'):
                profile_1d[j,i]   = hbt.nanmean(im[is_good])
        
            if (method == 'median'):
                profile_1d[j,i] = hbt.nanmedian(im[is_good])
    
            # clipped = astropy.stats.sigma_clip(im[is_good], sigma=1.5)
            
            # std_1d[j,i] = np.nanstd(clipped)
        
    if do_std: 
        return (radius_1d, profile_1d, std_1d)
    
    return (radius_1d, profile_1d)

### END OF FUNCTION DEFINTION
    
if (__name__ == '__main__'):
    pass



### END OF FUNCTION DEFINTION
    
if (__name__ == '__main__'):
    pass
