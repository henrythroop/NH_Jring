#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 15:24:45 2018

@author: throop
"""

def plot_img_wcs(img, wcs, ra=274.73238190163573, dec=-20.86594759658748, title=None):

    # Default RA, Dec above is for MU69 on approach
    
    import matplotlib.pyplot as plt # pyplot
    import astropy.visualization

    """
    Just plots an image, and uses WCS to plot a specific RA/Dec on it.
    
    One-off func; not in a class.
    """

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    plt.imshow(stretch(img))
    (x, y) = wcs.wcs_world2pix(ra, dec, 0)
    plt.plot(x, y, marker = 'o', ms = 5, color='red')
    plt.title(title)
    plt.show()
    