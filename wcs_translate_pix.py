#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 16:21:43 2018

@author: throop
"""
import numpy as np
import astropy
from   astropy.wcs import WCS

from astropy.io import fits
import hbt
from plot_img_wcs import plot_img_wcs




def wcs_translate_pix(wcs, dx_pix, dy_pix):
    """
    Function takes a WCS header, and offsets it by a specified number of pixels.
    
    This changes the value of CRPIX (ie, which pixel is defined as the center).
    At first I was modidfying only CRVAL (ie, the RA/Dec val of the center point), but that 
    had was only 'mostly' accurate, and applying a negative translation didn't totally undo it.
    
    The WCS passed is modified in place. There is no return value.
    
    This function is standalone -- not part of any class.
    
    Note that if I translate by 100 pixels, and then back by -100 pixels, I do not reach 
    exactly the starting location. It is close, but not exact. I think this is OK, but I'm
    not totally sure the reason.
     
    'x' is horizontal direction (as plotted by matplotlib). Positive means to translate to the right.
    'y' is vertical   direction.                            Positive means to translate up (with origin=0 set).
    
    Parameters
    -----
    
    dx_pix, dy_pix:
        Amount to translate by, in pixels in the x and y directions. x and y are vertical and horizontal, 
        as per matplotlib images.
    
    """

    print(f'Translating WCS coords to offset by dx={dx_pix}, dy={dy_pix}')
    wcs.wcs.crpix[0] += dx_pix
    wcs.wcs.crpix[1] += dy_pix

def wcs_zoom(wcs, zoom, shape_orig):
    """
    Zoom the WCS by a certain factor.
    e.g., zoom=4 → make each pixel into 4x4=16.
    
    The WCS passed is modified in place. There is no return value.
    
    This function is standalone -- not part of any class.
    
    This changes the value of CRPIX (ie, the pixel value of the center point), and pc + cd (the pixel scales)
    
    """
    
    try:
        wcs.wcs.pc = wcs.wcs.pc / zoom
    except AttributeError:
        wcs.wcs.cd = wcs.wcs.cd / zoom   # Four-element pixel scale matrix
    
    crpix = wcs.wcs.crpix  # Location x and y of center. x = horizontal.
    
    # Now set the pixel location of the center, stored in 

    # XXX not sure if this is quite correct yet!
    
    crpix_new = np.array( [(crpix[0]+0.5)/shape_orig[0] * (zoom * shape_orig[0]) - 0.5,
                           (crpix[0]+0.5)/shape_orig[1] * (zoom * shape_orig[1]) - 0.5] )
    
    wcs.wcs.crpix = crpix_new
    
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
    hdulist = fits.open(file)
    arr = hdulist['PRIMARY'].data
    hdulist.close()
    wcs = WCS(file)

    hbt.figsize((12,12))    
    plot_img_wcs(arr,wcs, title='Original')
    
    dx_pix = 12
    dy_pix = 27
    
    print(f'Center = {wcs.wcs.crval}')
    wcs_translate_pix(wcs, dx_pix, dy_pix)
    print(f'Center = {wcs.wcs.crval}')

    print()
    
    wcs_translate_pix(wcs, -dx_pix, -dy_pix)
    print(f'Center = {wcs.wcs.crval}')
    
    # Now, read in an image, pad it, translate it, and see if I can keep the WCS properly done.

    wcs = WCS(file)
    
    plot_img_wcs(arr,wcs, title='Original')

    dxy_pad = 10
    arr_pad = np.pad(arr, ((dxy_pad, dxy_pad),(dxy_pad, dxy_pad)), 'constant')
    wcs_pad = wcs.deepcopy()
    wcs_translate_pix(wcs_pad, dxy_pad, dxy_pad)
    plot_img_wcs(arr_pad, wcs_pad, title=f'Padded by {dxy_pad}')
    
    arr_pad_trans = np.roll(np.roll(arr_pad, dx_pix, axis=1), dy_pix, axis=0)
    wcs_pad_trans = wcs_pad.deepcopy()
    wcs_translate_pix(wcs_pad_trans, dx_pix, dy_pix)
    plot_img_wcs(arr_pad_trans, wcs_pad_trans, title=f'Padded by {dxy_pad}')
    
    