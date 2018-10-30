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
    
    This changes the value of CRVAL (ie, the RA/Dec val of the center point). Nothing else is changed.
    
    The WCS passed is modified in place. There is no return value.
    
    This function is standalone -- not part of any class.
    
    Note that if I translate by 100 pixels, and then back by -100 pixels, I do not reach 
    exactly the starting location. It is close, but not exact. I think this is OK, but I'm
    not totally sure the reason.
     
    'x' is horizontal direction (as plotted by matplotlib).
    'y' is vertical   direction
    
    Parameters
    -----
    
    dx_pix, dy_pix:
        Amount to translate by, in pixels in the x and y directions. x and y are vertical and horizontal, 
        as per matplotlib images.
    
    """

    print(f'Translating WCS coords to offset by dx={dx_pix}, dy={dy_pix}')
    crpix = wcs.wcs.crpix  # Center pixels. Often 127.5, 127.5
    crval = wcs.wcs.crval  # Center RA, Dec. In degrees
    
    wcs_original = wcs.deepcopy()
    
    # Calculate the new center position (in RA/Dec, based on the pixel shifts)
    
    radec_shifted = wcs.wcs_pix2world(np.array([[crpix[0] + dx_pix, crpix[1] + dy_pix]]), 0)
    
    # Set the new center position, in radec.
    
    wcs.wcs.crval = np.ndarray.flatten(radec_shifted)
    
    wcs_new = wcs.deepcopy()
    
    # Now do a check. Get RA/Dec of center in old, and check pixel position of this RA/Dec in new
    
    radec_original = wcs_original.wcs_pix2world([crpix],0) # ← use 0 for pix2world
    pix_shifted = wcs_new.wcs_world2pix(radec_original,1)  # use 1 for world2pix. See documentation -- not clear.
    
    print(f'In original, center at RA/Dec {crval}')
    print(f'In original, center at Pixel  {crpix}')
    print(f'In new, this RA/Dec is at pixel location {pix_shifted}')
    
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
    
# =============================================================================
# Define a test function which is called when we run this file. This is just an example of using the function.
# =============================================================================
    
if (__name__ == '__main__'):
    
    file = '/Users/throop/Data/MU69_Approach/porter/KALR_MU69_OpNav_L4_2018301/lor_0403016638_0x633_pwcs2.fits'
    
    wcs = WCS(file)
    
    dx_pix = 00
    dy_pix = 50
    
    print(f'Center = {wcs.wcs.crval}')
    wcs_translate_pix(wcs, dx_pix, dy_pix)
    print(f'Center = {wcs.wcs.crval}')

    print()
    
    wcs_translate_pix(wcs, -dx_pix, -dy_pix)
    print(f'Center = {wcs.wcs.crval}')
    
