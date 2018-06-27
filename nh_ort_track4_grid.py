#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 00:05:08 2018

@author: throop
"""


import pdb
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
import os.path
import os

import astropy
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library

import pickle # For load/save

import gzip

from   datetime import datetime

# HBT imports

import hbt

# =============================================================================
# Define a class to handle all of the details of the Track-4 grids.
# This class does not create the grids, but it does read, write, plot, and fly thru.
# The ultimate goal of this class is to create the time-dependent particle densities that
# get given to Doug Mehoke for Track 5.
#
# THE CURRENT FILE ONLY HAS CLASS DEFINITIONS FOR THE TRACK4 GRIDS.
# IT HAS NO EXECUTABLE CODE.
#
# To run Track 4:
#
#   - Execute nh_track4_calibrate.py . This reads in all of DPH/DK's individual dust trajectories,
#     and merges them into '4D' dust grids, which are properly calibrated to match a given I/F.
#     Typically this reads in 108 files, and outputs 64 files, named *.grids4d.gz. These grids
#     are essentially just matrices (7, 200, 200, 200) with the dust density as a func of XYZ and grain size.
#
#   - Then execute nh_ort_track4_flyby.py. This reads all of the 64 grids files, and 
#     outputs a list of dust densities vs. time, for each one. 
#     Output is a table, essentially showing dust density (in # km-3) as a func of grain size, and time.
#     Typically 64 files, *.dust .
#
# The class definition and function are in the current file (NH_ORT_TRACK4_GRID.PY).
#     
# =============================================================================
    
class nh_ort_track4_grid:
    
    def __init__(self, arg):
    
        """
        Save the grid to the class, and do any initialization.
        
        Takes one input. It can be either of:
            
            density: A 4D grid: (num_particle_sizes, x, y, z). Typically (7, 200, 200, 200).
            
        or
          
            file: A filename, in `.grid4d.gz` format.
        
        """
        
        if (type(arg) == str):
            self.read(arg)

        else:
            self.density = arg
            
        self.num_grids = hbt.sizex(self.density)
        self.name_trajectory = None                 # Set as default, but will be changed later     
        self.name_test       = 'ort4-ring'         
        self.axis_sum = 1   # Which axis do we sum along for images? Should be as visible from Sun and/or SC.
                            # 1 → Sum along Y dir, which is Sun dir. This will make sunflower rings visible.  
            
        
# =============================================================================
# Set all ring parameters based on passed-in values
# =============================================================================
    
    def set_parameters(self, 
                             albedo=None, speed=None, q=None, rho=None,   # These are the physical params,
                                                                          # and they go into the filename and DPH sims.
                             b=None, beta=None, s=None,                   # These are the particle sizes we sum over.
                                                                          # They are used for plotting, but not else. 
                             name_trajectory=None, 
                             name_test=None,
                             resolution_km=None):
    
        """
        Save the parameters corresponding to a run.

        I guess this is kind of a getter-setter function. It's really not necessary. But it does make 
        for one easy line of code to set all variables, I guess.
        
        parameters:
            albedo, speed, q, rho    -- These are for the runs. One per 4D grid.
            b, beta, s    -- These are for particle size, one per grid.
            
            *** Note that for Python logic to work, the  arrays b, beta, s must be passed as lists, 
                 and *not* NumPy arrays. Testing against None works on lists, but chokes on NumPy arrays.
        """
                                       
        if albedo:
            self.albedo = albedo

        if speed:
            self.speed = speed
                                                                                      
        if q:
            self.q = q

        if rho:
            self.rho = rho
            
        # Now process the keywords related to particle size. These should each have the same length as .num_grids.
        
        if b:
            self.b = b

        if beta:
            self.beta = beta
                                                                                      
        if s:
            self.s = s

        if name_test:
            self.beta = beta
                                                                                      
        if name_test:
            self.s = s
            
        if resolution_km:
            self.resolution_km = resolution_km

# =============================================================================
# Make plots of flattened grids, one plot for each particle size
# =============================================================================
                                                         
    def plot(self, stretch_percent=98, axis_sum = None):
        
        """
        Make a plot of all of the sub-grids in a grid. Each one is flattened, and they are all ganged up.
        
        Optional parameters
        ---
        
        axis_sum:
            0=X, 1=Y, 2=Z
            
        """

# Here is how the order of plotting works.
# - Grid start out as a 3D array, in XYZ.
# - If we flatten in one axis, then the other two are left. e.g., flatten in Y, and then array is an XZ array.
# - Then when we imshow(), the first axis is vertical, and second is horizontal. e.g., X on vert, Z on horizontal.
# - If we do a transpose() on the array before imshow(), then order above is reversed (first is on horizont, etc).
# - Finally, if we pass the 'origin=lower' keyword, then the vertical axes will be drawn bottom-up, not top-down.
        
        if axis_sum is None:
            axis_sum = self.axis_sum
            
        (_, nx, ny, nz) = np.shape(self.density)
        
        xpos = hbt.frange(-nx/2, nx/2)*self.resolution_km[0]
        ypos = hbt.frange(-ny/2, ny/2)*self.resolution_km[1]
        zpos = hbt.frange(-nz/2, nz/2)*self.resolution_km[2]

        extent = (ypos[0], ypos[-1], zpos[0], zpos[-1])
        
        origin = 'lower'  # Need this to flip vertical axis properly
        
        extent = np.array(extent)/1000   # Truncate to units of 1000 km, not km.
        
        stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
#        hbt.figsize((15,15))
        for i in range(self.num_grids):
            b_i     = np.amax(self.b) - i     # b goes in opposite order from size, so plot it backwards.
            beta_i  = self.beta[i]              # Calc beta so we can get size. MRS slide 5.5
            s_i     = self.s[i]                 # Particle size. MRS slide 6.4. Millimeters.
   
            plt.subplot(1, self.num_grids, i+1)
            plt.imshow(np.transpose(stretch(np.sum(self.density[i], axis=axis_sum))), extent=extent, origin=origin)
            plt.title(r's={:.2f} mm, $\beta$={:.1e}'.format(s_i, beta_i))
            plt.tight_layout()

            if (axis_sum == 0):  # If summing along X
                plt.ylabel('Z [Mm]')   # We *always* use np.transpose() above, so that reverses usual axis plot order!
                plt.xlabel('Y [Mm]')
            if (axis_sum == 1):  # If summing along Y
                plt.ylabel('Z [Mm]')
                plt.xlabel('X [Mm]')
            if (axis_sum == 2):  # If summing along Z
                plt.ylabel('Y [Mm]')
                plt.xlabel('X [Mm]')

        plt.show()
        print(f'  albedo={self.albedo}, q={self.q}, rho={self.rho}, speed={self.speed}')


# =============================================================================
# Make a plot showing the geometry of the path thru the system
# =============================================================================

    def plot_trajectory_geometry(self):
    
        """
        Make a plot showing the geometry of the path thru the system
        Six plots, showing Radius, Velocity, X Y Z position, etc.
        """
        
        hbt.figsize((8,5))
        plt.subplot(2,3,1)
        plt.plot(self.delta_et_t, np.array(self.radius_t)/1000)
        plt.ylim((0,np.amax(self.radius_t)/1000))
        plt.title('Radius')
        plt.xlabel('dt from CA [sec]')
        plt.ylabel('Radius [Mm]')
        
        plt.subplot(2,3,2)
        plt.plot(self.delta_et_t, np.array(self.lat_t) * hbt.r2d)
        plt.title('Lat [deg]')
        
        plt.subplot(2,3,3)
        plt.plot(self.delta_et_t, np.array(self.lon_t) * hbt.r2d)
        plt.title('Lon [deg]')
        
        plt.subplot(2,3,4)
        plt.plot(self.delta_et_t, np.array(self.x_t)/1000)
        plt.title('X')
        plt.xlabel('dt from CA [sec]')
        plt.ylabel('Radius [Mm]')
        
        plt.subplot(2,3,5)
        plt.plot(self.delta_et_t, np.array(self.y_t)/1000)
        plt.title('Y')
        
        plt.subplot(2,3,6)
        plt.plot(self.delta_et_t, np.array(self.z_t)/1000)
        plt.title('Z')
        
        plt.tight_layout()
                    
        plt.show()

# =============================================================================
#  Make a Q&D plot of line-of-sight optical depth through the grid
# =============================================================================
        
    def plot_tau(self, stretch_percent=98, axis_sum = 1):
        """
        Do a Q&D sum to plot line-of-sight optical depth through the grid. 
        
        Sum in the Y dir, which is the line-of-sight from sun or observer.
        
        Units of 'density' are particles per km3  [MRS slide 6.6] 
        
        Optional inputs
        ----
        
        stretch_percent:
            Amount to scale by
            
        axis_sum:
            Axis along which to sum. 0=X, 1=Y, 2=Z.
            
            
        """

        origin = 'lower'  # Make sure the plot is oriented the right way for imshow.
        
        (_, nx, ny, nz) = np.shape(self.density)
        
        xpos = hbt.frange(-nx/2, nx/2)*self.resolution_km[0]
        ypos = hbt.frange(-ny/2, ny/2)*self.resolution_km[1]
        zpos = hbt.frange(-nz/2, nz/2)*self.resolution_km[2]

        extent = np.array((ypos[0], ypos[-1], zpos[0], zpos[-1]))/1000

        self.calc_tau()
        
        do_stretch_linear = False

        img = self.tau_2d
        
        # Stretch the image, as appropriate
        
        val = np.percentile(img[img > 0], 50)  # Set a constant offset for the log stretch

        img_stretch = hbt.logstretch(img, val)
        if do_stretch_linear:
            img_stretch = astropy.visualization.PercentileInterval(98)(img)
      
        # Create a colorbar, in the plot itself
        
        img_stretch_t = np.transpose(img_stretch)  # Transpose the image, swapping X and Y axes, so Z is vertical axis.
        
        colorbar = hbt.frange(np.amin(img_stretch_t), np.amax(img_stretch_t), hbt.sizey(img_stretch_t))
        img_stretch_t[:,-1] = colorbar  # There is probably a better way to do this?
        img_stretch_t[:,-2] = colorbar
        img_stretch_t[:,-3] = colorbar
        img_stretch_t[:,-4] = colorbar
        img_stretch_t[:,-5] = colorbar

        # Render the image
        
        plt.imshow(img_stretch_t, origin=origin, extent=extent)
                
        # Create the labels for the colorbar, and individually place them
            
        num_ticks_colorbar = 5 # Number of vertical value to put on our colorbar
        
        for j in range(num_ticks_colorbar):
            range_stretch = hbt.frange(np.amin(img_stretch_t), np.amax(img_stretch), num_ticks_colorbar)
            val_text = hbt.logstretch_invert(range_stretch[j], val)
#            val = round(val)  # Convert to zero if it is very close
            xval = 0.55*np.max(extent)
            yrange = np.max(extent)-np.min(extent)
            yval = np.min(extent) + 0.02*yrange + (j/(num_ticks_colorbar-1) * 0.92*yrange)
            if not(do_stretch_linear):
                plt.text(xval, yval, f'{val_text:.2e}', color = 'white') # Label the colorbar, if log stretch only.
                    
        plt.title(f'Tau, Merged, {len(self.s)} sizes')
        
#        tau_(max,typ) = ({self.tau_max:.1e}, {self.tau_typ:.1e}), ' + 
#                  f'I/F_(max,typ) = ({self.iof_max:.1e}, {self.iof_typ:.1e})')
        
        plt.xlabel('X [1000 km]')
        plt.ylabel('Z [1000 km]')
        plt.show()




# =============================================================================
# Calculate optical depth tau and iof in various ways through the grid
# =============================================================================
        
    def calc_tau(self):
        """
        Calculate optical depth tau and iof in various ways through the grid
        
        Note that this tau (and I/F) is summed only over a limited set of particle size bins.
        So, this tau will necessarily be less than the tau of the ring that we 
        orignally were trying to match. That is just how MRS's eq's 6.4 + 6.6 work.
        
        Save results as variables in the method:
            tau_2d
            tau_max
            iof_max
            tau_typ
            iof_typ
            
        """
        
        tau = np.zeros((200,200))

        for j,s in enumerate(self.s):  # Loop over particle size 's' in mm
                        
            tau_i = np.sum(self.density[j], axis=1) * (math.pi * s**2) * 1e-12  # units of mm2/km2, so then convert
            tau += tau_i
#            print(f'Summing with size = {s} mm')
        
        tau *= 250   # Adjust for fact that bin density is # per km^3, but bin size is 250 km
        
        tau_max     = np.amax(tau)
        iof_max     = tau_max * self.albedo
        tau_typ     = np.percentile(tau, 99)
        iof_typ     = tau_typ * self.albedo
        
        # Save values to the class, if we need them later
        
        self.iof_max = iof_max
        self.iof_typ = iof_typ
        self.tau_max = tau_max
        self.tau_typ = tau_typ
        
        self.tau_2d = tau
        
# =============================================================================
# Create the output filename automatically, based on the parameter values
# This is the filename used for Doug Mehoke.       
# =============================================================================

    def create_filename_track4(self, base = 'ort4'):
        
        str_traj = self.name_trajectory
        
        str_base = base
        
        str_speed = 'y{:3.1f}'.format(abs(self.speed))
        
        str_q = 'q{:3.1f}'.format(abs(self.q))
        
        str_albedo = 'pv{:4.2f}'.format(self.albedo)
        
        str_rho = 'rho{:4.2f}'.format(self.rho)
        
        file_out = f"{str_base}_{str_traj}_{str_speed}_{str_q}_{str_albedo}_{str_rho}.dust"
        
        return file_out
    
# =============================================================================
# Write out the 4D grid to disk, in a pickle format that can be used to fly through the ring.
# My code will later read this disk, and use the output of that to give to Doug Mehoke.        
# =============================================================================
    
    def write(self, file=None, dir=None, do_compress=False, style = 'pickle'):
        """
        Write the 4D grid itself to a file. The run parameters (albedo, rho, q, speed) are encoded into the filename.
        The remaining parameters (b, beta, s) are written.
        
        This file is a temporary file which we will then fly through, using nh_ort_track4_flyby().
        
        Format of the file is a pickle tuple, with fields:
              (density, albedo, rho, q, speed, b, beta, s, (resolution_km)).
              
        'density' is the large 4D array, with format (num_grids, x, y, z) -- typically (7, 200, 200, 200).      
        
        ** Warning: These are really big files! 400+ MB by default. Will 'compress' down to 7 MB.
        
        https://stackoverflow.com/questions/18474791/decreasing-the-size-of-cpickle-objects
                
        Optional Parameters:
            file: Filename. By default, it is auto-generated from the parameters.
            
            dir:  Directory
            
            do_compress: 
                Boolean. Do we use gzip compression on the output grids, or not?
                
        """
        
        if not dir:
            dir = os.path.expanduser('~/Data/ORT4/throop/track4/')

        if not file:
            file = self.create_filename_track4()
        
    # If requested, write the output to a pickle file (optionally compressed),
    # which we use for the next phase of the hazard proc, when we fly thru it for Doug Mehoke.
    
        if (style == 'pickle'):
            
            file = file + '.pkl'
            
            if do_compress:
                
                file = file + '.gz'
                
                protocol=-1    
                with gzip.open(os.path.join(dir,file), 'wb') as f:
                    pickle.dump((self.density, 
                                 self.albedo, self.rho, self.q, self.speed, 
                                 self.b, self.beta, self.s, 
                                 self.resolution_km), 
                                 f, 
                                 protocol)
    
            else:
                lun = open(os.path.join(dir, file), 'wb')        
                pickle.dump((self.density, 
                                 self.albedo, self.rho, self.q, self.speed, 
                                 self.b, self.beta, self.s, 
                                 self.resolution_km), lun)
                lun.close()
     
            print("Wrote: " + os.path.join(dir,file))
 
    # If requested, write the output to a .xyz file, for use in MeshLab.
    # This is more of a 'sparse' file, since it simply lists the indices of each filled grid.
    # No densities, no grids with zeros, etc.
    # If a grid has >= 1 particle  in it, it's listed.
    # If it     has    0 particles it it, it's omitted.

        if (style == 'xyz'):
            
            file = file + '.txt'  # Extension .xyz sometimes works, but .txt is def more reliable
            
            # Collapse the density array, and sum over particle size
            
            density = np.sum(self.density,axis=0)
            (nx, ny, nz) = np.shape(density)            
            
            # Take a threshhold -- that is, only list the the points above a density %ile of __.
            # If we list all the points, it is simply too dense a grid.
                        
            thresh_low  = np.percentile(density[density > 0], 90)
            thresh_high = np.percentile(density[density > 0], 97)

            str_red     = '255; 0; 0'
            str_green   = '0; 255; 0'
            str_blue    = '0; 0; 255'
            str_yellow  = '255; 255; 0'
            
            f = open(os.path.join(dir,file), 'w')

            (num_pts_high, num_pts_low) = (0,0)

            # Write the file itself, by looping over every grid site, and printing a line if value is large enough
            
            for i_x in range(nx):
                for i_y in range(ny):
                    for i_z in range(nz):

                          if thresh_low < density[i_x, i_y, i_z] < thresh_high:
                              num_pts_low += 1
                              f.write(f'{i_x/nx - 0.5:0.2f}; {i_y/ny - 0.5:0.2f}; {i_z/nz - 0.5:0.2f}; {str_red}\n')

                          if density[i_x, i_y, i_z] > thresh_high:
                              num_pts_high += 1
                              f.write(f'{i_x/nx - 0.5:0.2f}; {i_y/ny - 0.5:0.2f}; {i_z/nz - 0.5:0.2f}; {str_yellow}\n')
                              
            f.write('1.00')  # I don't know what this is. But the file seems to require it.
            f.close()
            print(f'Wrote: {os.path.join(dir,file)}, {num_pts_low}+{num_pts_high} points')
                
        return

# =============================================================================
# Read a 4D grid file from disk
# =============================================================================
    
    def read(self, file, dir=None):
        """
        Read a 4D grid file from disk

        Format of the file is a pickle tuple, with (grid_4d, albedo, rho, q, v, b, beta, s).
        
        File may be either .grid4d or .grid4d.gz . Reader will decompress if needed.
        
        Optional Parameters:
            dir: Directory
            
        """
        
        if not dir:
            dir = os.path.expanduser('~/Data/ORT4/throop/track4/')

        print("Reading: " + file)

        if 'gz' in file:
            with gzip.open(os.path.join(dir, file), 'rb') as f:
                 (self.density, 
                                 self.albedo, self.rho, self.q, self.speed, 
                                 self.b, self.beta, self.s, 
                                 self.resolution_km) = pickle.load(f)
            f.close()
        else:
            f = open(os.path.join(dir, file), 'rb')
            (self.density, 
                             self.albedo, self.rho, self.q, self.speed, 
                             self.b, self.beta, self.s, 
                             self.resolution_km) = pickle.load(f)            
            f.close() 
        return        
    
# =============================================================================
# Fly a trajectory through the grids and sample it
# =============================================================================
        
    def fly_trajectory(self, name_observer, et_start, et_end, dt):

        """
        Sample the ring along a flight path.
        
        Parameters
        ----
        
        et_start:
            
        et_end:
        
        dt:
            
        """
    
        # Save the passed-in parameters in case we need them 
        
        self.et_start = et_start
        self.et_end   = et_end
        self.dt       = dt
        
        self.name_observer = name_observer
        
        self.area_sc  = 5*u.m**2              # Area of spacecraft
        
        self.abcorr = 'LT'
        
        # Set up the output time array
        
        num = math.ceil( (et_end - et_start) / dt.to('s').value ) + 1
        
        et_t = hbt.frange(int(et_start), int(et_end), num)
        
        # Calc offset in DT from C/A time
            
        delta_et_t = et_t - np.mean(et_t)

        # Set up position arrays for center of bins
        
        self.x_1d = hbt.frange(-hbt.sizex(self.density[0])/2, hbt.sizex(self.density[0])/2-1)*self.resolution_km[0]
        self.y_1d = hbt.frange(-hbt.sizey(self.density[0])/2, hbt.sizey(self.density[0])/2-1)*self.resolution_km[0]
        self.z_1d = hbt.frange(-hbt.sizez(self.density[0])/2, hbt.sizez(self.density[0])/2-1)*self.resolution_km[0]
        
# Set up output arrays for all of the quantities. These will be time series.
        
        radius_t = []
        x_t      = []
        y_t      = []
        z_t      = []
        v_t      = [] # Velocity (scalar)
        lon_t    = []
        lat_t    = []
        density_t= []
        bin_x_t  = []
        bin_y_t  = []
        bin_z_t  = []
                
        # Loop over et

        for i,et_i in enumerate(et_t):

    # Get the position of s/c in MU69 frame. We just look up XYZ position, in the frame.
    # Then we figure out what box that is in, and look up the density there.
    # This is really very simple, with no PXFORM required. All of it happens automatically after we 
    # set the frame to be the sunflower frame.
              
            (st, lt) = sp.spkezr(self.name_observer, et_i, self.frame, self.abcorr, self.name_target)
                                                                    # Get position of s/c in MU69 frame!
                                                                    
            (radius, lon, lat) = sp.reclat(st[0:3])
            
            # Get the lon/lat wrt time. Note that for MU69 flyby, the lat is always positive.
            # This is non-intuituve, but it is because the MU69 Sunflower frame is defined s.t. the ring
            # is in the XZ plane. 
            # In this plane, indeed NH stays 'above' MU69 the whole flyby, with lat always positive.
            
            radius_t.append(radius)
            lon_t.append(lon)
            lat_t.append(lat)

            ### If requested, send the s/c right down the midplane, at impact parameter zero.
            ### This is just for testing, and should never be used in production.
            
            DO_OVERRIDE_XZ = False
            
            if (DO_OVERRIDE_XZ):
                st[0] = 0
                st[2] = 0
                
                print("Warning: Setting X = Z = 0 for trajectory. Do not use this in production!!")
                
            # Get the XYZ positions wrt time.
        
            x_t.append(st[0])
            y_t.append(st[1])
            z_t.append(st[2])
            
            v_t.append( sp.vnorm(st[3:6]) )
        
        # End of loop over t. Now that we have done the flythru, we look up all the index values.
        
        # Convert into bin values. This is vectorized.
        # ** If values exceed the largest bin, then the index returned will be too large for the density lookup!
        
        bin_x_t = np.digitize(x_t, self.x_1d)
        bin_y_t = np.digitize(y_t, self.y_1d)
        bin_z_t = np.digitize(z_t, self.z_1d)
        
        v_t = np.array(v_t) * u.km/u.s
        
        # If indices are too large, drop them by one. This just handles the edge cases so the code doesn't crash.
        
        bin_x_t[bin_x_t >= len(self.x_1d)] = len(self.x_1d)-1
        bin_y_t[bin_y_t >= len(self.y_1d)] = len(self.y_1d)-1
        bin_z_t[bin_z_t >= len(self.z_1d)] = len(self.z_1d)-1
        
        # Now that we have all the XYZ bins that the s/c travels in, get the density for each one.
        # We should be able to do this vectorized -- but can't figure it out, so using a loop.
        
        number_t = np.zeros((self.num_grids, len(et_t)))
        
        for i in range(len(et_t)):
            number_t[:,i] = self.density[:,bin_x_t[i], bin_y_t[i], bin_z_t[i]]  # Number per km3

        # Calculate cumulative number sum. Just for fun.
        
        number_cum_t       = np.cumsum(number_t, axis=1)
    
        # Now for fun, calculate the number of grains that intercept a s/c of a given area.
        # We do this calc very crudely, just by taking the area of s/c vs area of a bin edge. We are ignoring
        # the slant path of the s/c, so we could underestimate the flux by up to sqrt(3).
    
        # Get the volume swept up by the sc at each timestep. It is in km3.
        
        vol_t = (self.area_sc * v_t * self.dt).to('km^3')
        
        # Get total number swept up: (# per km3) (volume of path)
        
        number_sc_t = number_t * vol_t
        
        number_sc_cum_t = np.cumsum(number_sc_t, axis=1)
            
        # Save all variables so they can be retrieved
        
        self.radius_t      = radius_t  # This is distance from the center, in 3D coords. Not 'ring radius.'
        self.et_t          = et_t
        self.delta_et_t    = delta_et_t
        self.bin_x_t       = bin_x_t
        self.bin_y_t       = bin_y_t
        self.bin_z_t       = bin_z_t
        self.number_t      = number_t
        self.number_sc_t   = number_sc_t
        self.number_sc_cum_t=number_sc_cum_t
        self.number_cum_t  = number_cum_t
        self.lon_t         = lon_t
        self.lat_t         = lat_t
        self.x_t           = x_t
        self.y_t           = y_t
        self.z_t           = z_t
        self.v_t           = v_t
        
        return
    
# =============================================================================
# Output the trajectory and particle intercept info to a file
# =============================================================================
        
    def output_trajectory(self, suffix=None, do_positions=True):

        """
        Write an output table. The filename is auto-generated based on the parameters already supplied.
        
        The output value is in particles per km3.
                            
        Parameters
        -----
        do_positions:
            If True, output three columns which describe the XYZ position of the spacecraft as function of time.
        str:
            A string to be inserted into the filename, as part of the filename.
        """

        # Choose the method of table writing. Astropy is easier to code, but I can't use my own header to match MRS's.
        
        do_write_astropy = False
        
        # Create the table
        
        # We really don't care about fractional values on any of these. So, truncate everything
        # to be an integer. There is no easy way to force all columns to print as int, so just do it manually.
        
        t = Table([np.array(self.et_t).astype(int), 
                   np.array(self.delta_et_t).astype(int)], 
                      names = ['ET', 'ET from CA'])

        if do_positions:
            t['X [km]'] = np.array(self.x_t).astype(int)
            t['Y [km]'] = np.array(self.y_t).astype(int)
            t['Z [km]'] = np.array(self.z_t).astype(int)
        
        for i in range(len(self.s)):    # Loop over particle size, and add a new table column for each particle size.
                t['n_{}, {:.3f}, # per km3'.format(i, self.s[i])] = np.array(self.number_t[i]).astype(int)
        
        # Create the output filename
        
        dir_out = '/Users/throop/Data/ORT4/throop/track4/'
        
        file_out = self.create_filename_track4()
        
#        file_out = 'ort1_q{}_i{}_a{}.txt'.format(self.q_dust, self.inclination, self.albedo)
#        if do_positions:
#            file_out = file_out.replace('.txt', '_pos.txt')
#       
        path_out = os.path.join(dir_out, file_out)
        
        if do_write_astropy:

            # Write the file using Astropy
            
            t.write(path_out, format = 'ascii', overwrite=True)
            print("Wrote: {} using astropy writer".format(path_out))

        else:

            # Write the file using explicit manual formatting
            # There are seven size radial bins.
            
            lun = open(path_out, "w")
            lun.write("#    First line is '0' and then size bins, in mm\n")
            lun.write("#    Remaining are delta_ET, n(r_1), n(r_2), n(r_3), n(r_4), n(r_5), n(r_6), n(r_7)\n")
            lun.write("#    n(r) are number per km3\n")
            lun.write("#    Henry Throop {}\n".format(str(datetime.now())))          
            lun.write("{} {} {} {} {} {} {} {}\n".format(0, 
                                                    self.s[0],
                                                    self.s[1],
                                                    self.s[2],
                                                    self.s[3],
                                                    self.s[4],
                                                    self.s[5],
                                                    self.s[6]))
            for i in range(len(t)):                                                
                lun.write("{} {} {} {} {} {} {} {}\n".format(
                        t[i][1],
                        t[i][2],
                        t[i][3],
                        t[i][4],
                        t[i][5],
                        t[i][6],
                        t[i][7],
                        t[i][8]))
            lun.close()
            print("Wrote: {} using manual writer".format(path_out))
                   
        return
    