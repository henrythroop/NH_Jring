#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 12:45:37 2017

Extracts radial *or* azimuthal profile from unwrapped ring images.
This is a more flexible routine than the one that extracts *both* profiles

  nh_jring_extract_profile_from unwrapped:  This routine,  more flexible. One profile at a time.
  nh_jring_extract_profiles_from unwrapped: Other routine, less flexible. Two profiles at a time.

bins_radius:  Input: the definition of radius  bins in the unwrapped image (required)
bins_azimuth: Input: The definition of azimuth bins in the unwrapped image (required)

range_profile [if azimuthal]:  Portion of radius  range to use, when extracting azimuthal profile
                       Can be a scalar fraction (0.5), or vector fraction (0.3, 0.6), or vector distance (128000,129000)
   
range_profile [if radial]:     Portion of azimuth range to use, when extracting radial    profile.
                       Can be a scalar fraction (0.5) only.

type_profile:        String, 'radial' or 'azimuthal' (case-insensitive, first letter checked only)

mask_unwrapped: A boolean mask, indicating whether to use individual pixels in the final output, or not.

@author: throop
"""

def nh_jring_extract_profile_from_unwrapped(im_unwrapped, 
                                            bins_radius, 
                                            bins_azimuth,
                                            range_profile, 
                                            type_profile,    # String: either 'radial' or 'azimuthal'
                                            mask_unwrapped=False):

    """ Extract a radial or azimuthal profile, as specified.
    
        :im_unwrapped:   The unwrapped image, N x M pixels
        :bins_radius:    Array defining the radial bins. M elements long.
        :bins_azimuth:   Array defining the azimuthal bins. N elements long.
        :range_profile:  Range of azimuth to use (for radial profile), or v/v
        :type_profile:   String: either 'radial' or 'azimuthal'. Case insensitive.
        :mask_unwrapped: 
            
            
    Parameters
    ----------
    im_unwrapped : 
        The unwrapped image, N x M pixels
    bins_radius :
        Array defining the radial bins. M elements long.
    third : 
        Array defining the azimuthal bins. N elements long.
    range_profile :
        Range of azimuth to use (for radial profile), or v/v
    type_profile :
        String: either 'radial' or 'azimuthal'. Case insensitive. 
    mask_unwrapped :
        
    Returns
    -------
    Array:
        The profile itself, either radial or azimuthal. The length is as specified in the input parameters.    
    
    """
        
    import hbt
    import numpy as np
    from   scipy.interpolate import griddata
    import warnings

    if type(mask_unwrapped) == type(np.array([])):
        DO_MASK = True

    if (DO_MASK):    
        is_good_unwrapped = (mask_unwrapped == True)     # Keep True values as it

        im_unwrapped_masked = im_unwrapped.copy()
        im_unwrapped_masked[is_good_unwrapped == False] = np.nan   # Convert False values into NaN in the image array
    
        im2 = im_unwrapped_masked

    else:
        im2 = im_unwrapped

# Extract the radial profile, if requested

    if (type_profile.upper()[0] == 'RADIAL'[0]):
        
# Extract the radial profile, if requested.

#  We use a subset of the azimuth angles.
#    e.g., if range=0.1, then extract the central 10% of the azimuth angles
#    e.g., if range=(0.9, -0.4), then extract the central 90% of azimuth angles, but then *exclude* 
#         the central 40% of angles          

# Take the profile itself. Use mean() along one axis. 
# We do not want to sum -- we want to use mean to preserve flux.
         
        # In case of a single range -- e.g., 0.1
        
        if isinstance(range_profile, float) or isinstance(range_profile, int):
            bin_0 = int(hbt.sizex(bins_azimuth) * (0.5 - (0.5 * range_profile)))  
            bin_1 = int(hbt.sizex(bins_azimuth) * (0.5 + (0.5 * range_profile)))
    
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)    
                profile_radius  = np.nanmean(im2[:, bin_0:bin_1],1)          # Use nanmean() which ignores NaN
                                                                             # Problem happening here        
        
        # In case of a double range -- e.g., (0.9, -0.3)
        
        elif isinstance(range_profile, tuple):
            bin_0 = int(hbt.sizex(bins_azimuth) * (0.5 - (0.5 * range_profile[0])))  
            bin_1 = int(hbt.sizex(bins_azimuth) * (0.5 - (0.5 * np.abs(range_profile[1])))) 
            bin_2 = int(hbt.sizex(bins_azimuth) * (0.5 + (0.5 * np.abs(range_profile[1])))) 
            bin_3 = int(hbt.sizex(bins_azimuth) * (0.5 + (0.5 * range_profile[0])))

            print("Extracting from bins {}:{}, {}:{}".format(bin_0, bin_1, bin_2, bin_3))
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)    
                profile_radius  = (np.nanmean(im2[:, bin_0:bin_1],1) +
                                   np.nanmean(im2[:, bin_2:bin_3],1)) / 2

        else:
            
            raise(ValueError("Cannot parse range! range_profile = {}".format(range_profile)))
          
        return(profile_radius)
        
# Extract the azimuthal profile, if requested
   
    if (type_profile.upper()[0] == 'AZIMUTHAL'[0]):
   
      if (isinstance(range_profile, float) or isinstance(range_profile, int)):   
                                        # If passed a single radial range value, like 0.5, then use as fraction
          bin_0 = int(hbt.sizex(bins_radius) * (0.5 - (0.5 * range_profile)))
          bin_1 = int(hbt.sizex(bins_radius) * (0.5 + (0.5 * range_profile)))

                                        # If passed two radial range values, like (0.3, 0.4), then use as fractions
      elif ( isinstance(range_profile, tuple) and (range_profile[0] < 1) and (range_profile[1] < 1) ):
          bin_0 = int(hbt.sizex(bins_radius) * (range_profile[0]))
          bin_1 = int(hbt.sizex(bins_radius) * (range_profile[1]))
    
                                        # If passed two radial range values, like (128000,129000), then use as dist's.
      elif ( isinstance(range_profile, tuple) and (range_profile[0] > 1) and (range_profile[1] > 1) ):
          bin_0 = hbt.x2bin(range_profile[0], bins_radius)
          bin_1 = hbt.x2bin(range_profile[1], bins_radius)

      else:
          raise(ValueError("Cannot parse range!"))
          
      with warnings.catch_warnings():
         warnings.simplefilter("ignore", category=RuntimeWarning)    
         
         # Take the profile itself. Use mean() along one axis. 
         
         profile_azimuth = np.nanmean(im2[bin_0:bin_1, :],0)           
         
      return(profile_azimuth)