#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 16:21:43 2018

@author: throop
"""
import numpy as np
import astropy
from   astropy.wcs import WCS

def wcs_translate_pix(wcs, dx_pix, dy_pix):
    """
    Function takes a WCS header, and offsets it by a specified number of pixels.
    
    This changes the value of CRVAL. Nothing else is changed.
    
    The WCS passed is modified in place. There is no return value.
    
    This function is standalone -- not part of any class.
    
    Parameters
    -----
    
    dx_pix, dy_pix:
        Amount to translate by, in pixels in the x and y directions.     
    
    """

    
    print(f'Translating WCS coords to offset by dx={dx_pix}, dy={dy_pix}')
    crpix = wcs.wcs.crpix  # Center pixels. Often 127.5, 127.5
    crval = wcs.wcs.crval  # Center RA, Dec. In degrees
    
    # Calculate the new center position (in RA/Dec, based on the pixel shifts)
    
    radec_shifted = wcs.wcs_pix2world(np.array([[crpix[0] + dx_pix, crpix[1] + dy_pix]]), 0)
    
    # Set the new center position, in radec.
    
    wcs.wcs.crval = np.ndarray.flatten(radec_shifted)
    
    # Q: Does PC array need to change? PC = pixel coord transformation matrix.
    # http://docs.astropy.org/en/stable/api/astropy.wcs.Wcsprm.html
    # Might be deprecated. Well, no. has_pc(): 
    #   "PCi_ja is the recommended way to specify the linear transformation matrix."
    # http://www.stsci.edu/hst/HST_overview/documents/multidrizzle/ch44.html -- this looks like 
    #   the CD matrix (which is similar to PC) specifies rotation, but not offsets.
    # PC1_1 = CD1_1 if CDelt = 1, which it is. 
    # http://docs.astropy.org/en/stable/api/astropy.wcs.Wcsprm.html#astropy.wcs.Wcsprm.has_cd
    # So, I should be good to leave PC alone.

#    return wcs