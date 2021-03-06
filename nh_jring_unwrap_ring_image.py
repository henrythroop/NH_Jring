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
                               binsize_azimuth, 
                               planes, dx=0, dy=0, mask_stray=None, mask_objects=None):
        
    """
    Unwrap a 2D ring image from (RA, Dec) into (radius, azimuth).
    
    Parameters:
    -----    
    
    All of these parameters are mandatory.
    
    im:      
        Image array
        
    num_bins_radius:
        Number of radial bins. **If this is too large, then we get a 'stretching' artifact**.
        
    limits_radius:
        Extent of the inner and outer radius, in km.
        
    binsize_azimuth:
        size of azimuthal bins, in radians.
        
    planes:  
        The table of backplanes
        
    mask_stray:    
        2D array of pixel flags. True = good. 
        
    mask_objects:
        2D array of pixel flags. True = good. 
  
    dx, dy: 
        Pixel values to roll the image by -- that is, an additional offset to be added to nav info in WCS header. 
          Both the image and the mask are rolled by this amount. Integer.
  
    Output
    ------
    (im_unwrapped, bins_radius, bins_azimuth)
  
    """  
        
    import hbt
    import numpy as np
    from   scipy.interpolate import griddata
    import math
    import matplotlib.pyplot as plt

# Process input

    if (mask_objects is None):    # NB: We need to use 'is' here. '==' will do an element-by-element comparison.
        DO_MASK_OBJECTS = False   # Alternatively, could also do 'if type(mask) is np.ndarray'
  
    else:
        DO_MASK_OBJECTS = True

    if (mask_stray is None):    # NB: We need to use 'is' here. '==' will do an element-by-element comparison.
        DO_MASK_STRAY = False   # Alternatively, could also do 'if type(mask) is np.ndarray'
  
    else:
        DO_MASK_STRAY = True
    
        
    # https://stackoverflow.com/questions/36783921/check-if-variable-is-none-or-numpy-array-in-python
    # https://stackoverflow.com/questions/15008380/double-equals-vs-is-in-python
    
#    DO_MASK = (mask != None)  # If we are passed the 'mask' value, then 
        
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
    
    # Quantize these values to the next-lower and next-upper bins (e.g., 0.001543 → 0.001), to use a common base.
    
    azimuth_seg_start = azimuth_seg_start - np.mod(azimuth_seg_start, binsize_azimuth) - binsize_azimuth
    azimuth_seg_end   = azimuth_seg_end   - np.mod(azimuth_seg_end,   binsize_azimuth) + binsize_azimuth

    # Calculate number of bins
    
    # NB: Must use round() not int() here to avoid 15.999999 → 15. And then int for benefit of numpy.
    
    num_bins_azimuth = int( round( (azimuth_seg_end - azimuth_seg_start) / binsize_azimuth) + 1 )
        
    indices_3_good = (azimuth_all_3 >= azimuth_seg_start) & (azimuth_all_3 < azimuth_seg_end)
    
    azimuth_all_good = azimuth_all_3[indices_3_good]
    radius_all_good  = radius_all_3[indices_3_good]
    dn_all_good      = dn_all_3[indices_3_good]
    
    # Extract arrays with the proper pixel values, and proper azimuthal values
    
    azimuth_all = azimuth_all_good
    radius_all  = radius_all_good
    dn_all      = dn_all_good
    
    if (DO_MASK_OBJECTS):
        mask_objects_roll = np.roll(np.roll(mask_objects, int(round(dy)), 0), int(round(dx)), 1)
        mask_objects_all = mask_objects_roll[is_ring_all]
        mask_objects_all_3 = np.concatenate((mask_objects_all, mask_objects_all, mask_objects_all))
        mask_objects_all_good    = mask_objects_all_3[indices_3_good]
        mask_objects_all    = mask_objects_all_good
        
    if (DO_MASK_STRAY):

        mask_stray_roll = np.roll(np.roll(mask_stray, int(round(dy)), 0), int(round(dx)), 1)     
        mask_stray_all = mask_stray_roll[is_ring_all]
        mask_stray_all_3 = np.concatenate((mask_stray_all, mask_stray_all, mask_stray_all))
        mask_stray_all_good    = mask_stray_all_3[indices_3_good]
        mask_stray_all    = mask_stray_all_good
    
#==============================================================================
#  Now regrid the data from xy position, to an unrolled map in (azimuth, radius)
#==============================================================================

# Construct the gridded image line-by-line

    dn_grid         = np.zeros((num_bins_radius, num_bins_azimuth))  # Row, column
    
    bins_azimuth    = hbt.frange(azimuth_seg_start, azimuth_seg_end, num_bins_azimuth)
    bins_radius     = hbt.frange(limits_radius[0], limits_radius[1], num_bins_radius)        

    if (DO_MASK_OBJECTS):
        mask_objects_grid     = dn_grid.copy() 
    
    if (DO_MASK_STRAY):
        mask_stray_grid       = dn_grid.copy() 
    
    for i in range(num_bins_radius-1):  # Loop over radius -- inner to outer. Do one radial output bin at a time.
        
        # Select only bins with right radius and azimuth
        
        is_ring_i = np.array(radius_all > bins_radius[i]) & np.array(radius_all < bins_radius[i+1]) & \
                    np.array(azimuth_all > azimuth_seg_start) & np.array(azimuth_all < azimuth_seg_end) 
        
        if np.sum(is_ring_i) > 0:
            dn_i         = dn_all[is_ring_i]  # Get the DN values from the image (adjusted by nav pos error)
            radius_i     = radius_all[is_ring_i]
            azimuth_i    = azimuth_all[is_ring_i]
            
            
            # Now make sure there is not a large gap in azimuth in the data. If there is, then make it out w/ NaN's.
            # This happens when the ring ansae gets near the image edge. griddata() will blindly interpolate,
            # but we want to fill gap w/ NaN so that it will not do so.

            do_plug_gaps = True
            
            if do_plug_gaps:
                
                # Put the azimuth points in monotonic order. They are not otherwise (and griddata does not require).
                # But, to look for gaps, we need to sort
               
                d_azimuth_gap_max = 0.05     # Max width of a gap, in radians. Larger than this will be NaN-plugged.
                
                indices   = np.argsort(azimuth_i)
                azimuth_i_sort = azimuth_i[indices]
                dn_i_sort      = dn_i[indices]
                
                # Check the point-by-point increase in azimuth
                
                d_azimuth_i_sort = azimuth_i_sort - np.roll(azimuth_i_sort,1)   
                
                # Flag it if it exceeds some pre-set size
                
                is_gap = np.abs(d_azimuth_i_sort) > d_azimuth_gap_max
                
                # And then put NaNs at the start and end of the gap
                
                is_gap[0] = False
                if np.sum(np.where(is_gap)):
                    w = np.where(is_gap)[0][0]
                    dn_i_sort[w-1:w+2] = np.nan

                grid_lin_i   = griddata(azimuth_i_sort, dn_i_sort, bins_azimuth, method='linear')
                dn_grid[i,:] = grid_lin_i         # Write a row of the output array   
            
            # Now do the interpolation. If there are NaNs it will skip that region.
            else:
                grid_lin_i   = griddata(azimuth_i, dn_i, bins_azimuth, method='linear')
                dn_grid[i,:] = grid_lin_i         # Write a row of the output array   

            if DO_MASK_STRAY:
                mask_stray_i         = mask_stray_all[is_ring_i]  # Get the DN values from the image 
                                                                  # (adj by nav pos error)
                grid_lin_i           = griddata(azimuth_i, mask_stray_i, bins_azimuth, method='linear')
                mask_stray_grid[i,:] = grid_lin_i

            if DO_MASK_OBJECTS:
                mask_objects_i         = mask_objects_all[is_ring_i]
                grid_lin_i             = griddata(azimuth_i, mask_objects_i, bins_azimuth, method='linear')
                mask_objects_grid[i,:] = grid_lin_i

# The way I've done the regridding, each output bin goes from n:n+1. So, the final bin is empty. Just duplicate it.
# There is probably a better way to do this.
# This kind of means that my output array is shifted 1/2-bin from my input array. Ugh.
                
# Also, this routine is quite slow, and I really think it introduces some artifacts. If I was a CS major, 
# I could write a much better version of this algorithm.               
#
# Algorithm uses griddata(), which is just interpolation. So, while it does not preserve 'area*flux', it does preserve
# the actual DN values ('DN per pixel') -- and it creates new pixels, with the same DN-per-pixel values.
# So, a radial profile, or azimuthal profile, by taking the mean along one axis, should work as intended, and give
# a value in DN-per-pixel.
                
    n                        = num_bins_radius
    dn_grid[n-1,:]           = dn_grid[n-2,:]
    mask_stray_grid[n-1,:]   = mask_stray_grid[n-2,:]
    mask_objects_grid[n-1,:] = mask_objects_grid[n-2,:]                    
                
# Save the variables, and return

    image_unwrapped    = dn_grid     # NB: This has a lot of NaNs in it.
    
    if (DO_MASK_OBJECTS and DO_MASK_STRAY):
        mask_stray_unwrapped   = mask_stray_grid
        mask_objects_unwrapped = mask_objects_grid
        
        # Convert these to booleans, and set any NaN to be False (since they are probably off-edge)
        
        
        return (image_unwrapped, mask_stray_unwrapped, mask_objects_unwrapped, bins_radius, bins_azimuth)
    
    else:
        return (image_unwrapped, bins_radius, bins_azimuth)
    

# =============================================================================
# Do some tests and validation of the unwrapping functions.
# These are good to have in general, but I wrote them to track down a specific error in J-ring navigation.
# Specifically, we want to unwrap the backplane, to see if we get what we should.
# =============================================================================

#%%
       
def test():

#%%
    
    import hbt
    import pickle
    import astropy
    import matplotlib.pyplot as plt
    import numpy as np
    from nh_jring_unwrap_ring_image import nh_jring_unwrap_ring_image
    from skimage.io import imread, imsave   # For PNG reading
#    from pyds9 import DS9

    file_pickle = '/Users/throop/Data/NH_Jring/out/nh_jring_read_params_571.pkl' # Filename to read to get files, etc.

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    lun = open(file_pickle, 'rb')
    t = pickle.load(lun)
    lun.close()

    # Process the group names. Some of this is duplicated logic -- depends on how we want to use it.

    groups = astropy.table.unique(t, keys=(['Desc']))['Desc']
    
    index_group = 8
#    index_image = [32,33,34,35]
#    index_image = [24,25]
    index_image = [99]

#    index_image = [32]
    file_mask = '/Users/throop/Data/NH_Jring/masks/mask_{}_{}.png'.format(index_group, 32)
    file_mask = '/Users/throop/Data/NH_Jring/masks/mask_{}_{}-{}.png'.format(index_group, 0,7)

    num_bins_radius = 2000  # NB: My results are indep of this. That is, E and W ansa are off by 600 km, regardless
                           # of the number of radial bins chosen. So it's not just a fencepost error.
    num_bins_azimuth = 5000
    
#    num_bins_radius = 300
#    num_bins_azimuth = 500
    
    limits_radius = np.array([120000, 130000])

    groupmask = (t['Desc'] == groups[index_group])
    t_group = t[groupmask]
    dir_backplanes = '/Users/throop/data/NH_Jring/out/'

    image = {}  # Set up a dictionary to store all my output images. "Keys of a dictionary may be strings or ints."
    image_raw = {} # Raw image, before any bg subtraction, or unwrapping, etc.
    planes = {} # Set up dictionary to store planes

    profile_radius = {}         # Radial profile of raw image
    profile_radius_masked = {}  # Profile of masked image 
    profile_azimuth = {}        # Azimuthal profile 
    profile_azimuth_masked = {}
    profile_radius_plane = {}
    profile_azimuth_plane = {}
    radius_unwrapped = {}
    image_unwrapped = {}
        
# Read the images into an array
    
    for i in index_image:
        
        # Read image from disk
        image_raw[i] = hbt.read_lorri(t_group[i]['Filename'], frac_clip = 1, bg_method = 'None')
        
        # Process it: remove sfit, load the stray light mask, etc.
        
        out = hbt.nh_jring_process_image(image_raw[i], 'String', 'p5 mask_7_0-7')
#        out = hbt.nh_jring_process_image(image_raw[i], 'String', '7/32-35 p5')

        # Grab the mask used during the sfit(), if it was given to us.
        
        if isinstance(out, tuple):
            image[i], mask_stray = out
        
        # Or otherwise, just get the image and no mask
        
        else:
            image[i] = out
     
            mask_stray = imread(file_mask) > 127   # Convert the PNG file to boolean
        
        file_backplane = dir_backplanes + t_group['Shortname'][i].replace('.fit', '_planes.pkl')
    
    # Save the shortname associated with the current backplane. 
    # That lets us verify if the backplane for current image is indeed loaded.

        file_short = t_group['Shortname'][i]
				
        lun = open(file_backplane, 'rb')
        print("Reading backplane {}".format(file_backplane))
        
        planes[i] = pickle.load(lun)
        lun.close()
        
# Now unwrap the image

        hbt.figsize((15,5)) 
        
        aspect = 1/30
        
        # XXX bug: the final row of every array returned by this function is all 0. Fix this.
        
        (image_unwrapped[i], mask_stray_unwrapped, mask_objects_unwrapped, bins_radius, bins_azimuth) \
                                = nh_jring_unwrap_ring_image(image[i], num_bins_radius, 
                                                     limits_radius, num_bins_azimuth,
                                                     planes[i], mask_stray = mask_stray,
                                                     mask_objects = mask_stray)
                                
        extent = [np.amin(bins_azimuth), np.amax(bins_azimuth), 
                  np.amax(bins_radius)/1000, np.amin(bins_radius)/1000]
    
        # Plot the unwrapped image
        
        hbt.figsize((25,25))
        plt.imshow(stretch(image_unwrapped[i]), extent=extent, aspect=aspect*3, cmap=plt.cm.plasma)
        plt.title("{}/{}: {}".format(index_group, i, file_short))
        plt.show()

        # Plot the unwrapped image, with mask overlaid
 
        plt.imshow(stretch(image_unwrapped[i]), extent=extent, aspect=aspect, cmap=plt.cm.Greys_r)

        plt.imshow(mask_stray_unwrapped, extent=extent, aspect=aspect, alpha=0.2, 
                   cmap=plt.cm.Reds_r,vmin=-0.5,vmax=0.5)

        plt.title("{}/{}: {} MASK OVERLAY".format(index_group, i, file_short))
        plt.show()        
        
        # Compute radial profile of the raw image
                
        profile_radius[i]  = np.nanmedian(image_unwrapped[i], axis=1)
        profile_radius[i] -= np.amin(profile_radius[i]) # Subtract off a DC signal from rad prof

        # Compute radial profile of the masked image (ie, with stray removed)
        
        # XXXX There is some error in here. It somehow relates to taking mean of nan's vs. 0's, that kind of thing.
        
        arr = image_unwrapped[i] * mask_stray_unwrapped
        arr[arr == 0] = np.nan                           # Convert any 0's to NaN. This is for edges, and masked.
        
        profile_radius_masked[i]  = np.nanmedian(arr, axis=1)
#        profile_radius_masked[i] -= np.amin(profile_radius_masked[i]) 
        
        # Plot radial profiles
        
        plt.plot(bins_radius, profile_radius[i])
        plt.plot(bins_radius, profile_radius_masked[i])
        
        # Unwrap the backplane itself, which tests the SPICE navigation
        
        DO_UNWRAP_BACKPLANE = False
        
        if DO_UNWRAP_BACKPLANE:
            (radius_unwrapped[i], bins_radius, bins_azimuth) \
                                    = nh_jring_unwrap_ring_image(planes[i]['Radius_eq'], num_bins_radius, 
                                                         limits_radius, num_bins_azimuth,
                                                         planes[i])
                                    
            # Compute radial profile of the backplane (just for validation -- should be smooth horizontal gradient)
            
            profile_radius_plane[i]  = np.nanmean(radius_unwrapped[i], axis=1)
            plt.imshow(stretch(radius_unwrapped[i]), extent=extent, aspect=aspect)
            plt.title("{}/{}: {} -- BACKPLANE ONLY".format(index_group, i, file_short))
            plt.show()
    
# =============================================================================
#  All computations done -- plot the two radial profiles, and the two backplane profiles
# =============================================================================

# Plot radial profiles of data
   
    hbt.figsize((10,6)) 
     
    alpha = 0.5
    
    for i in index_image:    

        plt.plot(bins_radius/1000, profile_radius[i],        linestyle = 'solid', marker = 'None', 
                 label = "{}/{}: {}".format(index_group, i, t_group['Shortname'][i]), alpha=alpha)
        
        plt.plot(bins_radius/1000, profile_radius_masked[i], linestyle = 'solid', marker = 'None', 
                 label = "{}/{}: {}".format(index_group, i, t_group['Shortname'][i] + ' MASKED'), alpha=alpha)

        
#    plt.ylim((000,20))
        
    plt.ylabel('DN')
    plt.xlim(hbt.mm(bins_radius/1000))
    plt.legend()
    plt.xlabel('Radius [1000 km]')
    plt.title('Radial Profiles, Data. Radial Bins = {}'.format(num_bins_radius))
    plt.show()


# Plot radial profiles of backplane proper
    
    if DO_UNWRAP_BACKPLANE:
        for i in index_image:    
            plt.plot(bins_radius/1000, profile_radius_plane[i]/1000, linestyle = 'dotted', marker = '.', label = repr(i))
            
        plt.ylabel('Mean Backplane Radius [1000 km]')
        plt.xlabel('Navigated Radius [1000 km]')
        plt.xlim((126.5,127.5))
        plt.ylim((126.5,127.5))
        plt.legend()
        plt.title('Radial Profiles, Backplane')
        plt.axvline(x=127)
        plt.axhline(y=127)
        plt.show()

# Make a scatterplot of the backplanes. For each pixel in the unwrapped image, plot its Y pos and its value. 

        for i in index_image:
            a = np.meshgrid(bins_azimuth, bins_radius)
            mesh_radius = a[1]
            plt.imshow(mesh_radius)
            plt.imshow(stretch(radius_unwrapped[i]))
            plt.title('Radius, Backplane, {}/{}'.format(index_group, i))
            plt.show()
        
        offset = 0.02
        
        for i in index_image:
            plt.plot(mesh_radius.flatten()/1000, radius_unwrapped[i].flatten()/1000 + offset*(i==23), 
                     linestyle='none', marker = '.', 
                      ms=1, label='{}/{}: {} {}'.format(index_group, i, t_group['Shortname'][i],
                                     'plus offset' if (i==23) else ''))
            plt.ylim(limits_radius)
            plt.ylim((126.7,127.3))
            plt.xlim((126.7,127.3))
            plt.axvline(x=127)
            plt.axhline(y=127)
    
        plt.xlabel('Radius from bins')
        plt.ylabel('Radius from backplane pixels')
        plt.title('Radius of data in backplane, when unwrapped')
        plt.legend()            
        plt.show() 
         
# Plot the raw images, but mask the pixels between in a certain radial range, based on backplane.
# This is again for testing the backplane accuracy        
#
# *** This turned out to show my problem. In the J-ring GUI when I overlay the ring, it is in the same 
#     location in E and W ansae. But in the data here, it is *not* -- it's shifted by a few pixels.
#     This was due to using wrong abcorr when creating the backplanes. 4-Aug-2017.

    if DO_UNWRAP_BACKPLANE:
        hbt.figsize((10,10))
        for i in index_image:
            image_i = image[i].copy()
            radius_i  = planes[i]['Radius_eq']
            range_radius_mask = [128.9,129]
            is_ring = np.logical_and(radius_i > range_radius_mask[0]*1000, radius_i <= range_radius_mask[1]*1000)
            print("I")
            image_i[is_ring] = np.amax(is_ring)
            plt.imshow(stretch(image_i))
            plt.show()


if (__name__ == '__main__'):
    test()
