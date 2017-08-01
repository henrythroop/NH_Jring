# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 09:01:30 2016

@author: throop
"""

##########
# Unwrap the ring image, based on coordinates provided by the backplane
# Turns it from a plot in (RA, Dec) to one in (Az, Radius)
#
# This routine does *not* preserve flux.
##########


def nh_jring_unwrap_ring_image(im, 
                               num_bins_radius, limits_radius, 
                               num_bins_azimuth, 
                               planes, dx=0, dy=0, mask = False):
        
    """
  im: image array
  radius: 1D array of radius, in km
  azimuth: 1D array of azimuth (radians)
  planes: the table of backplanes
  
  dx, dy: Pixel values to roll the image by. Also rolls the mask.
  
  output: (im_unwrapped, radius, azimuth)
  
    """  
        
    import hbt
    import numpy as np
    from scipy.interpolate import griddata
    import math

# Process input

    if type(mask) == type(np.array([0])):  # If we are passed the 'mask' value, then 
        DO_MASK = True
        
# Extract fields from the backplane array
    
    radius  = planes['Radius_eq']           # These are properly 4x4 or 1x1, as desired.
    azimuth = planes['Longitude_eq']
    phase   = planes['Phase']
    
#==============================================================================
# Examine backplane to figure out azimuthal limits of the ring image
#==============================================================================

    bins_radius = hbt.frange(limits_radius[0], limits_radius[1], num_bins_radius)
        
    # Select the ring points -- that is, everything inside the mask
    
    is_ring_all = ( np.array(radius > limits_radius[0]) & np.array(radius < limits_radius[1]))
    
# Make sure there are some valid data points. If not, return an error
    
    if (np.sum(is_ring_all) == 0):
        print("Error: No valid ring points between radius {:,.0f} .. {:,.0f} km".format(limits_radius[0], 
                                                                                        limits_radius[1]))
        print("Image radius limits = {:,.0f} .. {:,.0f} km".format(np.amin(radius), np.amax(radius)))
        raise ValueError('NoValidRingPoints')
        
    radius_all  = radius[is_ring_all]     # Make a list of all of the radius values
    azimuth_all = azimuth[is_ring_all]  # Make a list of all of the azimuth points for all pixels
    phase_all   = phase[is_ring_all]

    dn_all     = np.roll(np.roll(im, int(round(dy)), 0), int(round(dx)), 1)[is_ring_all]   
                                                                # Axis 0 is y. 1 is x.
                                                                # If a float is passed, round to the closest int
        
#    phase_mean  = np.mean(phase_all)     # Get the mean phase angle across the ring.
#    
    # Now take these raw data, and rearrange them so that we can take the longest continuous segment
    # We do this by appending the timeseries to itself, looking for the largest gap (of no az data), 
    # and then the data will start immediately after that.
    
    # _2 indicates a double-length array (ie, with [azimuth, azimuth + 2pi])
    # _s indicates sorted
    # _d indicates delta
    
    azimuth_all_3 = np.concatenate((azimuth_all, azimuth_all + 2*math.pi, azimuth_all + 4*math.pi))
    dn_all_3      = np.concatenate((dn_all,      dn_all,                  dn_all))
    radius_all_3  = np.concatenate((radius_all, radius_all,               radius_all))
    
    azimuth_all_3_s = np.sort(azimuth_all_3, kind = 'heapsort')
    azimuth_all_3_s_d = azimuth_all_3_s - np.roll(azimuth_all_3_s, 1)
    
    # Look for the indices where the largest gaps (in azimuth) start
    
    # XXX we get an error here on images if there are no valid ring points.
    
    index_seg_start_3_s = (np.where(azimuth_all_3_s_d > 0.999* np.max(azimuth_all_3_s_d)))[0][0]
    index_seg_end_3_s   = (np.where(azimuth_all_3_s_d > 0.999* np.max(azimuth_all_3_s_d)))[0][1]-1
    
    # Get proper azimithal limits. We want them to be a single clump of monotonic points.
    # Initial point is in [0, 2pi) and values increase from there.
                                                       
    azimuth_seg_start = azimuth_all_3_s[index_seg_start_3_s] # Azimuth value at the segment start
    azimuth_seg_end   = azimuth_all_3_s[index_seg_end_3_s]   # Azimuth value at the segment end
    
    indices_3_good = (azimuth_all_3 >= azimuth_seg_start) & (azimuth_all_3 < azimuth_seg_end)
    
    azimuth_all_good = azimuth_all_3[indices_3_good]
    radius_all_good  = radius_all_3[indices_3_good]
    dn_all_good      = dn_all_3[indices_3_good]
    
    # Extract arrays with the proper pixel values, and proper azimuthal values
    
    azimuth_all = azimuth_all_good
    radius_all  = radius_all_good
    dn_all      = dn_all_good

    if (DO_MASK):
        mask_all = mask[is_ring_all]
        mask_all_3 = np.concatenate((mask_all, mask_all, mask_all))
        mask_all_good    = mask_all_3[indices_3_good]
        mask_all    = mask_all_good
        
#==============================================================================
#  Now regrid the data from xy position, to an unrolled map in (azimuth, radius)
#==============================================================================

# Construct the gridded image line-by-line

    dn_grid         = np.zeros((num_bins_radius, num_bins_azimuth))  # Row, column
    bins_azimuth    = hbt.frange(azimuth_seg_start, azimuth_seg_end, num_bins_azimuth)
    bins_radius     = hbt.frange(limits_radius[0], limits_radius[1], num_bins_radius)        

    if (DO_MASK):
        mask_grid       = dn_grid.copy() 
    
    for i in range(num_bins_radius-1):  # Loop over radius -- inner to outer. Do one radial output bin at a time.
        
        # Select only bins with right radius and azimuth
        
        is_ring_i = np.array(radius_all > bins_radius[i]) & np.array(radius_all < bins_radius[i+1]) & \
                    np.array(azimuth_all > azimuth_seg_start) & np.array(azimuth_all < azimuth_seg_end) 
        
        if np.sum(is_ring_i) > 0:
            dn_i         = dn_all[is_ring_i]  # Get the DN values from the image (adjusted by nav pos error)
            radius_i     = radius_all[is_ring_i]
            azimuth_i    = azimuth_all[is_ring_i]
            grid_lin_i   = griddata(azimuth_i, dn_i, bins_azimuth, method='linear')
            dn_grid[i,:] = grid_lin_i         # Write a row of the output array   

            if DO_MASK:
                mask_i       = mask_all[is_ring_i]  # Get the DN values from the image (adjusted by nav pos error)
                grid_lin_i   = griddata(azimuth_i, mask_i, bins_azimuth, method='linear')
                mask_grid[i,:] = grid_lin_i
                
# Save the variables, and return

    image_unwrapped    = dn_grid     # NB: This has a lot of NaNs in it.
    
    if (DO_MASK):
        mask_unwrapped = mask_grid
        return (image_unwrapped, mask_unwrapped, bins_radius, bins_azimuth)
    
    else:
        return (image_unwrapped, bins_radius, bins_azimuth)    