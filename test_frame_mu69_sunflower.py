#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 12:39:40 2018

@author: throop
"""

import spiceypy as sp
import hbt

# =============================================================================
# DEFINE TEST_FRAME_MU69_SUNFLOWER()
# =============================================================================

def test_frame_mu69_sunflower():
    
    """
    This is just a quick test routine to check that MU69 Sunflower frame is more-or-less working.
    
    This program prints the RA/Dec of MU69's Y and Z axes, under a variety of different conditions.
    
    No options.
    
    Things to verify here:
        - That the Z RA/Dec point toward the orbit pole RA/Dec in all cases. 
              Small deviations << 1 deg allowed
              
        - That the Y RA/Dec point toward the Sun always for _ROT Frames.
              (ie, Y RA/Dec changes with time.)
              
        - That the Y RA/Dec point toward the Sun for the _INERT frames on 1 Jan 2015,
              and slowly move after / before that, at roughly 1 deg/year.
              (ie, Y RA/Dec is pretty much fixed)
              
    16-Jan-2018. HBT verified that output values look pretty good.

    """
    
    tms = ['kernels_sunflower.tm']  # Define the metakernel, which in turn calls the sunflower .tf frame

    frames = ['2014_MU69_SUNFLOWER_ROT', '2014_MU69_SUNFLOWER_INERT']
        
    utcs = ['1 Jan 2005 00:00:00', '1 Jan 2015 00:00:00']
    
    frame_j2k = 'J2000'
    
    # Get values from Simon Porter. See email from MRS ~16-Jan-2018.
    # These RA/Dec values are also put into the .tf file which I have made.
    
    ra_mu69_sp = 272.426110231801*hbt.d2r
    dec_mu69_sp = 68.831520928192*hbt.d2r
    
    # Define the MU69 Z axis. We will rotate this, and it should point in specified direction.
    
    z_mu69 = [0, 0, 1] # MU69 +Z axis. It should point to the specified RA/Dec
    y_mu69 = [0, 1, 0] # MU69 +Y axis. It should point away from the Sun
    
    print("Simon Porter pole position:")
    print("RA = {}, Dec = {}".format(ra_mu69_sp*hbt.r2d, dec_mu69_sp*hbt.r2d))
    print("---")
    
    # Loop over input parameters. For each combination, do a calculation, and print the output.
    
    for tm in tms:    
        for frame in frames: 
            for utc in utcs:
                sp.furnsh(tm)
                et = sp.utc2et(utc)
                mx = sp.pxform(frame, frame_j2k, et)
            
                z_mu69_j2k = sp.mxv(mx, z_mu69)                
                (_, ra_z, dec_z) = sp.recrad(z_mu69_j2k)

                y_mu69_j2k = sp.mxv(mx, y_mu69)                
                (_, ra_y, dec_y) = sp.recrad(y_mu69_j2k)
                
                print("Metakernel: {}".format(tm))
                print("UTC:        {}".format(utc))
                print("Frame:      {}".format(frame))
                print("Matrix:     \n{}".format(mx))
                print("Y RA = {}, Dec = {}".format(ra_y*hbt.r2d, dec_y*hbt.r2d))
                print("Z RA = {}, Dec = {}".format(ra_z*hbt.r2d, dec_z*hbt.r2d))
                print("\n---\n")
        
# =============================================================================
# End of function definition
# =============================================================================
        
if (__name__ == '__main__'):
    test_frame_mu69_sunflower()