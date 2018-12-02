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
    
    Parameters
    -----
    
    wcs:
        The input WCS
    
    zoom: 
        The amount of zoom. 1 = no zoom. 2 = double in size, 1x1 → 2x2.
    
    shape_orig:
        np.shape() result from the original, non-zoomed array.
        
    """
    
    # First change the pixel scale
    
    try:
        wcs.wcs.pc = wcs.wcs.pc / zoom
    except AttributeError:
        wcs.wcs.cd = wcs.wcs.cd / zoom   # Four-element pixel scale matrix
    
    crpix = wcs.wcs.crpix  # Location x and y of center. x = horizontal.
    
    print(f'Before: crpix={crpix}')
    
    # Now change the center position
    # This math is a bit weird here! It gives what appears to be the right answer, but it has more terms in it 
    # than I'd expect. This answer is 'right' in that when I zoom an image with scipy.ndimage.zoom(), and then 
    # do this transformation here on the WCS, I get the same results. But, it is basically experimental.
    # 
    # The reason for all this is that WCS coords start, apparently, at 0 in python arrays, but 1 in FITS files,
    # or something bizarre like that. The documentation mentioned this. Here I am using both of course.
    # So maybe that is contributing to the confusion.
    #
    # Also, it might be that my arguments to the shift or the ndimage.zoom() are incorrect, and here I am 
    # doing weird manipulations to make up for my errors over there. That is basically OK though -- 
    # as long as I am consistent and get the same answers both ways in the end.

    crpix_new = np.array( [crpix[0] * zoom,
                          crpix[1] * zoom] )
    
    # crpix_new = np.array( [(crpix[0]+0.5)/shape_orig[0] * (zoom * shape_orig[0]) + 1 - (3/2.) * zoom,
                           # (crpix[1]+0.5)/shape_orig[1] * (zoom * shape_orig[1]) + 1 - (3/2.) * zoom] )
    
    print(f'After: crpix={crpix_new}')
    
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
    
    file = '/Users/throop/Data/MU69_Approach/throop/backplaned/' + \
           'KALR_MU69_OpNav_L4_2018335/lor_0405944938_0x633_pwcs2_backplaned.fits'
    hdulist = fits.open(file)
    arr = hdulist['PRIMARY'].data
    header = hdulist['PRIMARY'].header
    et =  header['SPCSCET']
    hdulist.close()
    wcs = WCS(file)
    hbt.figsize((12,12))    
    
    plot_img_wcs(arr,wcs, title='Original',
                 name_observer='New Horizons', name_target='MU69', et=et)
    
    # dx_pix = 12
    # dy_pix = 27
    
    # print(f'Center = {wcs.wcs.crval}')
    # wcs_translate_pix(wcs, dx_pix, dy_pix)
    # print(f'Center = {wcs.wcs.crval}')

    # print()
    
    # wcs_translate_pix(wcs, -dx_pix, -dy_pix)
    # print(f'Center = {wcs.wcs.crval}')
    
    # # Now, read in an image, pad it, translate it, and see if I can keep the WCS properly synced with the image.

    # wcs = WCS(file)
    
    # dxy_pad = 10
    # arr_pad = np.pad(arr, ((dxy_pad, dxy_pad),(dxy_pad, dxy_pad)), 'constant')
    # wcs_pad = wcs.deepcopy()
    # wcs_translate_pix(wcs_pad, dxy_pad, dxy_pad)
    # plot_img_wcs(arr_pad, wcs_pad, title=f'Padded by {dxy_pad}', 
    #              name_observer='New Horizons', name_target='MU69', et=et)

    # Now test zooming. Start from scratch
    
    # hdulist = fits.open(file)
    # arr_orig = hdulist['PRIMARY'].data
    # hdulist.close()
    # wcs = WCS(file)
    
    arr = arr_orig.copy()
    plot_img_wcs(arr, wcs, title='original',
                 name_observer='New Horizons', name_target='MU69', et=et)

    # Zoom

    zoom = 1    
    wcs_zoomed = wcs.deepcopy()
    wcs_zoom(wcs_zoomed, zoom, np.shape(arr_pad))
    arr_zoom = scipy.ndimage.zoom(arr, zoom)
    plot_img_wcs(arr_zoom, wcs_zoomed, title=f'Zoomed by {zoom}',name_observer='New Horizons', 
                 name_target='MU69', et=et, width=100)

    # Translate, then zoom
    zoom = 1
    dx_pix = 10
    dy_pix = 30
    wcs_zoomed = wcs.deepcopy()
    arr_trans = np.roll(np.roll(arr, dx_pix, axis=1), dy_pix, axis=0)
    wcs_zoom(wcs_zoomed, zoom, np.shape(arr_pad_trans))
    wcs_translate_pix(wcs_zoomed, dx_pix, dy_pix)
    plot_img_wcs(arr_trans, wcs_zoomed, title=f'Translate by {dx_pix}, {dy_pix}, zoomed by {zoom}',
                 name_target='MU69', et=et, width=100)

    # Translate, then zoom.

    zoom = 6
    dx_pix = -5
    dy_pix = 5
    wcs2 = wcs.deepcopy()
    arr2 = arr.copy()
    
    plot_img_wcs(arr_trans, wcs_zoomed, title=f'Translate by {dx_pix}, {dy_pix}, zoomed by {zoom}',
             name_observer='New Horizons', name_target='MU69', et=et, width=72)
    
    # Roll and zoom the image
    
    arr2 = np.roll(np.roll(arr2, dx_pix, axis=1), dy_pix, axis=0)  # axis=1, x, right.   axis=0, y, up
    arr2 = scipy.ndimage.zoom(arr2,zoom)
    
    # Translate and zoom the WCS
    
    wcs_translate_pix(wcs2, dx_pix, dy_pix)
    wcs_zoom(wcs2, zoom, np.shape(arr))
    plot_img_wcs(arr2, wcs2, title=f'Translate by {dx_pix}, {dy_pix}, zoomed by {zoom}',
             name_observer='New Horizons', name_target='MU69', et=et, width=800)
                 
    
    # Translate
    # This is a bunch of consecutive translations and zooms. It's doing pretty well.
    
    zoom = 0.5

    arr2 = np.roll(np.roll(arr2, dx_pix, axis=1), dy_pix, axis=0)
    arr2 = scipy.ndimage.zoom(arr2,zoom)
    wcs_translate_pix(wcs2, dx_pix, dy_pix)
    wcs_zoom(wcs2, zoom, np.shape(arr2))
    plot_img_wcs(arr2, wcs2, title=f'Translate by {dx_pix}, {dy_pix}, zoomed by {zoom}',
                 name_observer='New Horizons', name_target='MU69', et=et)

    dx_pix = -10
    dy_pix = -30        
    arr2 = np.roll(np.roll(arr2, dx_pix, axis=1), dy_pix, axis=0)
    wcs_translate_pix(wcs2, dx_pix, dy_pix)
    plot_img_wcs(arr2, wcs2, title=f'Translate by {dx_pix}, {dy_pix}, zoomed by {zoom}',
                 name_observer='New Horizons', name_target='MU69', et=et, width=100)
