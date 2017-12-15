#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:21:09 2017

Modified from calc_area_curve in salt_analyze_spectra.py

@author: throop
"""

import numpy as np

def area_between_line_curve(data_x, data_y, limits_x, bins=False, binning=None, mean=False):
    """ Calculate the area between a curve, and a straight line.
    
    The line is determined by specifiying the x positions of the start and end points.
    
    Value returned is typically positive for emission line, and negative for absorption.
    
    Parameters
    ---
    data_x:
        X array (e.g., wavelength)
    data_y:
        Y array (e.g., intensity)
    limits_x:
        Tuple or two-element array. Min and max values to use.
    bins:
        Boolean. If set, the values of LIMITS_X are integer bins, not x values.
        
    Output
    ----
    (area, line_out, bins)
    
    area:
        Scalar area. Positive if curve is above line.
    line_out:
        A fit to the line, covering the whole input data range
    bins:
        Two-element array, with the bin numbers of the starting and ending limits
   
    """
    
    # Look up the proper bin values, if requested
    
    if not bins:
        limits_bins = ( np.min(np.where(data_x > limits_x[0])), 
                        np.max(np.where(data_x < limits_x[1])) ) 
    else:
        limits_bins = limits_x
    
    bin_0 = limits_bins[0]
    bin_1 = limits_bins[1]    

    # Smooth the data, if requested. XXX not currently active
    
    data_y_smooth = data_y
    
    # Extract just the relevant region of the arrays
    
    data_x_use = data_x[bin_0:bin_1]
    data_y_use = data_y_smooth[bin_0:bin_1]

    dy = data_y_use[-1] - data_y_use[0]		# Get the overall dx and dy from one end of the band to the other
    dx = data_x_use[-1] - data_x_use[0]

    # Fit the line to the endpoints
    
    m = dy / dx		          # Slope
    b = data_y_use[0] - m * data_x_use[0] # Y-intercept

    line = m * data_x_use + b			# Get an equation for the line, so we can plot it.
    						# We don't fit to the line or anything.

    # Sum the curve
    
    area = sum(data_y_use - line)
    
    # Take the mean, if requested
    
    if mean:
        area = area / len(data_x)
        
    # Create the output line
    
    line_out = m * data_x + b    # A fit to the line, covering the whole input data range
    
    return (area, line_out, limits_bins)