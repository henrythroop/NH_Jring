#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 15:24:45 2018

@author: throop
"""

def plot_img_wcs(img, wcs, ra=274.73344, dec=-20.86170, markersize=5, alpha=0.5, 
                 color='red', title=None, do_zoom_center = False, do_show=True,
                 width=None, **kwargs):
 

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
        Width of plot, in pixels. It is assumed to be square (height=width). If missing, the entire plot is shown.
    
    do_show:
        If False, suppress the final plt.show(). This allows other planes to be plotted on top later if requested.
        
    """

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    plt.imshow(stretch(img), origin='lower', **kwargs)
    (x, y) = wcs.wcs_world2pix(ra, dec, 0)
    plt.plot(x, y, marker = 'o', markersize = markersize, alpha=alpha, color=color)

    # Set the width and height appropriately
    
    if not width:                        # If no width passed, then just use the default.
        width = np.shape(img)[0]
        
    hw = int(width/2)
    ycenter = int((np.shape(img)[0])/2)  # Vertical
    xcenter = int((np.shape(img)[0])/2)  # Horizontal
    
    plt.xlim((xcenter-hw, xcenter+hw))
    plt.ylim((ycenter-hw, ycenter+hw))
    
    # Make some markers on the edge, so we can tell if MU69 is actually properly centered, or not.
    
    plt.plot([ycenter, ycenter, ycenter+hw, ycenter-hw], [xcenter-hw, xcenter+hw, xcenter, xcenter], 
             marker = 'o', color='red', ls='None')
        
    plt.title(title)

    # Display the whole thing
    
    if do_show:
        plt.show()
    