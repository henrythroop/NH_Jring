#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 15:24:45 2018

@author: throop
"""

import spiceypy as sp
import astropy
from   astropy.io import fits
import numpy as np
import spiceypy as sp
from   astropy.visualization import wcsaxes
import hbt
from   astropy.wcs import WCS

def plot_img_wcs(img, wcs, ra=274.73344, dec=-20.86170, markersize=5, alpha=0.5, 
                 color='red', title=None, do_zoom_center = False, do_show=True,
                 width=None, name_observer='New Horizons', name_target = 'MU69', 
                 et = None, do_stretch=True,
                 do_inhibit_axes=False,
                 do_colorbar=False,
                 do_plot_fiducial=True,
                 label_axis=None,
                 # vmin=None, vmax=None,
                 scale_km=None,
                 **kwargs):
 

    # Default RA, Dec above is for MU69 on approach
    # (274.73344, -20.86170) = 2-Nov-2018. But it changes v slowly, so should be OK.
    
    import matplotlib.pyplot as plt # pyplot
    import astropy.visualization
    import numpy as np

    """
    Just plots an image, and uses WCS to plot a specific RA/Dec on it.
    By default, the RA/Dec is the inbound asymptote for NH MU69. 
    This RA/Dec should be good to within 1 pixel for 15-Aug .. 15-Dec 2018.
    
    This is a one-off func; not in a class.
    
    **kwargs are passed to imshow(), not to plot().
    
    Parameters
    -----
    
    do_zoom_center:
        Flag. Zoom in on just the center region, or show the entire plot.
        
    width:
        If set, crop the plot to the specified width in pixels, centered at the plot center.
        Plot is assumed to be square (height=width).
        If not set, the entire plot is shown.
    
    do_show:
        If False, suppress the final plt.show(). This allows other planes to be plotted on top later if requested.
    
    do_stretch:
        [Default True]. Boolean. If set, apply a stretch to the outuput image before plotting.
    
    do_inhibit_axes:
        Boolean. If set, do not print the x and y values on the respective axes.
    
    do_plot_fiducial:
        Boolean. If set, plot little red dots on the edges, showing the original image center location.
        
    ra: 
        RA to center on [degrees]
        
    dec:
        Declination to center on [degrees]
        
    name_target:
    name_observer:
    et:
        If *all* of these parameters are supplied, then the RA and Dec are taken from a call to SPKEZR,
        rather than from the supplied RA/Dec.
    
    scale_km:
        pixel scale in km/pix. If this is passed, then xy axes are labeled in km from center.
        If this is 'None', then xy axes are labeled in pixels.

    label_axis:
        Text to label the axes with (e.g., 'pixels' or 'km from MU69'). If None, then no label is plotted.
        
    """

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.


    
    # Set the width and height appropriately
    
    if not width:                        # If no width passed, then just use the default.
        width = np.shape(img)[0]
        
    hw = int(width/2)
    ycenter = int((np.shape(img)[0])/2)  # Vertical
    xcenter = int((np.shape(img)[0])/2)  # Horizontal
    
    x0 = xcenter-hw  # center positions, in pixels
    x1 = xcenter+hw
    y0 = ycenter-hw
    y1 = ycenter+hw

    # Crop the image if needed

    img_crop = img[x0:x1, y0:y1]

    ycenter_crop = int((np.shape(img_crop)[0])/2)  # Vertical
    xcenter_crop = int((np.shape(img_crop)[0])/2)  # Horizontal
    
    # Stretch the image vertically. Must do this after the crop, not before!
    
    vmax = np.percentile(img_crop, 95)
    vmin = np.percentile(img_crop,  5)    

    img_crop_stretch = np.clip(img_crop, a_min=vmin, a_max=vmax)

    # If requested, calculate the x and y axes, in km
    
    extent = None # Set ehe default value
    
    size = np.shape(img_crop_stretch)

    if scale_km is not None:
        extent = [-size[0]/2 * scale_km, size[0]/2 * scale_km,
                  -size[1]/2 * scale_km, size[1]/2 * scale_km]    
        # print(f'extent = {extent}')
        
    else:
       extent = [0, size[0], 0, size[1] ]

    # Do the actual imshow()
    
    if not(do_stretch):    
        fig = plt.imshow(img_crop, origin='lower', extent=extent, **kwargs)
    else:
        fig = plt.imshow(img_crop_stretch, origin='lower', extent=extent, **kwargs)
    
    # If requested, inhibit printing values on the x and y axes
    
    if do_inhibit_axes:
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        
    # If requested, use SPICE to determine the position to plot
    
    if name_target and name_observer and et:
        abcorr = 'LT'
        frame = 'J2000'
        (st, lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
        vec = st[0:3]
        (_, ra, dec) = sp.recrad(vec)
        ra *= hbt.r2d
        dec *= hbt.r2d
        
    (x, y) = wcs.wcs_world2pix(ra, dec, 0)    # The '0' is the 'origin': 1 for FITS, where [1,1] is origin.
                                              #                          0 for numpy, where [0,0] is origin.
                                              # (Though if I go into DS9, it says that LL pixel is 0,0.)
                                              # I assume that I should use 0 here

# Definition of 'origin', which is third argument of wcs_world2pix:  
#  "Here, origin is the coordinate in the upper left corner of the image. In FITS and Fortran standards, 
#   this is 1. In Numpy and C standards this is 0." 

    # print(f'Position xy: {x}, {y}')
    
    
    # plt.xlim((xcenter-hw, xcenter+hw))
    # plt.ylim((ycenter-hw, ycenter+hw))
    
    # Make some markers on the edge, so we can tell if MU69 is actually properly centered, or not.
    
    if do_plot_fiducial:
        if scale_km is not None:
            # Mark center
            # plt.plot(x, y, marker = 'o', markersize = markersize, alpha=alpha, color=color)  
            # Mark four edges
            plt.plot([0*scale_km, 0*scale_km, hw*scale_km, -hw*scale_km], [-hw*scale_km, hw*scale_km, 0,0], 
                 marker = 'o', color='red', ls='None')
        else:
            # Mark center
            # Mark four edges
            # plt.plot(0, 0, marker = 'o', markersize = markersize, alpha=alpha, color=color)              
            plt.plot([ycenter_crop, ycenter_crop, ycenter_crop+hw, ycenter_crop-hw], 
                     [xcenter_crop-hw, xcenter_crop+hw, xcenter_crop, xcenter_crop], 
                 marker = 'o', color='red', ls='None')

    # Label axes, if requested
    if label_axis is not None:
        # if scale_km is not None:        
            # text = 'km'
        # else:
            # text = 'pix'    
        plt.xlabel(label_axis)
        plt.ylabel(label_axis)
            
    plt.title(title)

    # Plot a colorbar, if requested
    
    if do_colorbar:
        plt.colorbar(fig, format='%.0e')
        
    # Display the whole thing
    
    if do_show:
        plt.show()

# =============================================================================
# Run the function as a test case, if requested.
# =============================================================================
        
if (__name__ == '__main__'):
    
    # Load an example FITS file
    
    file = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_Hazard_L4_2018325/' + \
           'lor_0405122518_0x633_pwcs2_backplaned.fits'
    
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
    
    # Read and process the file.
    
    hdu = fits.open(file)
    et  = hdu['PRIMARY'].header['SPCSCET']
    w   = WCS(file)
    img = hdu['PRIMARY'].data
    hdu.close()
    
    # Pass the image data, WCS, and ET, and plot it
    
    plot_img_wcs(img, w, name_observer='New Horizons', name_target = 'MU69', et=et, do_colorbar=True)

    plot_img_wcs(img, w, name_observer='New Horizons', name_target = 'MU69', et=et, scale_km = 100, 
                 do_colorbar=True, label_axis='km', cmap='plasma')



    plot_img_wcs(img_superstack_median_iof, wcs_superstack, cmap=cmap_superstack, 
#                 title = f'{str_stack}, {str_reqid}, ring at {a_ring_km} km', 
                 title = f'{str_stack}, {str_reqid}', 
                 # width=width_pix_plot,
                 do_show=False,
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                 # extent = None,
                 scale_km = 100,
                 label_axis='km',
                 do_colorbar=True, do_stretch=True)

    
#%%%    
    # d = hbt.dist_center(100)
    # plt.imshow(d[0:50,0:50],extent=[-1,1,-2,2], origin='lower', aspect=0.2)
    # # plt.xlim((0,50))
    # # plt.ylim((0,50))
    # # plt.plot([1,2],[3,5])
    # plt.show()