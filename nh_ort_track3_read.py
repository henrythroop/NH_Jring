#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:38:48 2018

This reads the NH MU69 ORT 'Track 3' data, which is the output of N-Body simulations done by
Doug Hamilton + David Kaufmann.

This code does not combine the runs, or make output files, or anything. 
This file *only* works with one individual run at a time. It does:
    
    - Reads the main output file
    - Read the header, which contains timestep info, # of grains, etc.
    - Write a trajectory to .xyz files which can be used to visualize in MeshLab.
    - Prints information about the run, from header
    - Plotting the trajectory from an individual file (in nicely project X, Y, Z plots)
    
Incorporates code fragment from DK, which does the actual file I/O.
 
@author: throop
"""

import math
import os.path
import os

import astropy
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
import numpy as np
import astropy.modeling
from   astropy import units as u           # Units library
import struct

# HBT imports

import hbt

class nh_ort_track3_read:

# =============================================================================
# This is the class to read, plot, and output the N-body 'grid files' produced by DPH and DK
# for NH KEM Hazard ORT's.
# =============================================================================
    
# =============================================================================
# Initialize the class
# Read the specified track3 file into memory, from disk
# NB: Large beta → Small grains        
# =============================================================================
    
    def __init__(self, name, dir_base='/Users/throop/data/ORT4/hamilton/deliveries', binfmt = 2):
        
        """
        Initialize the object. Read in data from a specified file, and do initial processing on it.
        
        Parameters
        ------
            
        name:
            The name of the run. This is a slash-separated string, such as 
            `'ort1-0001/1000010/speed1/beta_2.6e-04/subset00'`. These slashes are not the best way 
            to do this, and it's not very compact. But the names are inconsistent, and this is the 
            most reliable way to deal with it for now.
        
        Optional parameters
        ------
        
        binfmt:
            Type of file to read. 1 = old file (pure grid). 2 = new file (sparse matrix). For ORT1 and beyond, 
            the default is `binfmt=2`.
            
        dir_base:
            The top-level directory that I download files from ixion into. Typically this ends in `deliveries`.
        """    
        
        do_kaufmann  = ('kaufmann' in dir_base)
        do_hamilton  = ('hamilton' in dir_base) 
        
        self.name = name
    
        parts = name.split('/')

        dir = name
        
        # Construct the filename for the 'header', which is a short text file. It is basically the same
        # as the entire directory tree, but with '/' → '_'.
      
        # Also, construct the filename for the data array itself. It is trivial.
        # As per the wiki, the data and header files should be identically named. They are not!
        # https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Integrations-to-Density
        
        if do_kaufmann:
            file_header = 'model.header'
            file_data   = f'model.array{binfmt}'

        if do_hamilton:
            file_header = 'header.txt'
            file_data   = f'grid.array{binfmt}'
                
        file_header_full = os.path.join(dir_base, dir, file_header)
        file_data_full   = os.path.join(dir_base, dir, file_data)

        self.file_header_full = file_header_full
        self.file_data_full   = file_data_full
        
        # Read the datafiles. Using code fragment from DK email 13-Feb-2018.

        #### START OF DK CODE FRAGMENT ####
        
        # Read header file
        
        f = open(file_header_full, 'r')
        text = f.read()
        f.close()
        lines = text.split('\n')
        
        self.header = lines
        
        nx = int(lines[0].split()[0])
        ny = int(lines[1].split()[0])
        nz = int(lines[2].split()[0])
        km_per_cell_x = float(lines[3].split()[0])
        km_per_cell_y = float(lines[4].split()[0])
        km_per_cell_z = float(lines[5].split()[0])
        
        # Read and process binary file
        
        density = np.zeros((nx, ny, nz), np.float32)
        f = open(file_data_full, 'rb')
        data = f.read()
        f.close()
                
        if binfmt == 1:
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        start = 4*(i*ny*nz + j*nz + k)
                        (density[k,j,i],) = struct.unpack('f', data[start:start+4])
        else: # binfmt == 2
            statinfo = os.stat(file_data_full)
            nentries = int(statinfo.st_size / 10) # record size = 10 bytes (3 2-byte ints + 4-byte float)
                                                  # int / int = float, in Python -- not sure why
                   
            for i in range(nentries):
                # Loop over list of non-empty cells. For each one, read xyz position (6 bytes), and density (4 bytes)
                start = 10*i
                (ix,iy,iz) = struct.unpack('hhh', data[start:start+6])   # h = short int
#                (density[iz-1,iy-1,ix-1],) = struct.unpack('f', data[start+6:start+10])  # Original order ZYX
                (density[ix-1,iy-1,iz-1],) = struct.unpack('f', data[start+6:start+10])   # HBT modified order XYZ, 
                                                                                          # to match wiki API.
                
    
        #### END OF DK CODE FRAGMENT ####
        
        self.density = density
        
        # Read a bunch of values from the header, and save them into the object
        
        self.km_per_cell_x = km_per_cell_x
        self.km_per_cell_y = km_per_cell_y
        self.km_per_cell_z = km_per_cell_z
        
        self.body_id   = self.get_header_val('body_id')
        self.duration  = self.get_header_val('duration')*u.year
        self.time_step = self.get_header_val('time_step')*u.s
        self.beta      = self.get_header_val('beta')
        self.speed     = self.get_header_val('speed')    # -3 or -2.2
        self.subset    = self.get_header_val('subset')
        self.grains    = self.get_header_val('grains')
        self.state_file= self.get_header_val('state_file')
        self.obj_file  = self.get_header_val('obj_file')
        
        return  # Return the object so we can chain calls (but, we don't really do that here)

# =============================================================================
# Search through the header, and find the value corresponding to a given keyword
# =============================================================================

    def get_header_val(self, keyword):
        
        """ 
        Get a header value from a list of header lines.
        
        As per https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Integrations-to-Density,
        the header is supposed to be a fixed format, with one specific value per line. However, I am building
        this function more robustly, in case the specs change in the future.
        
        Parameters
        -----
        
        keyword:
            A string that matches the keyword for a specified entry in the header. `keyword` should be chosen so that
            it will not match multiple lines in the header.
        
        """
       
                                                 
        out = ''
        for line in self.header:                    # Loop over every line
            if keyword in line:
                out = line.split('#')[0].strip()
                if hbt.is_number(out):                  # Convert from string to float, if possible
                    out = float(out)
        return out

# =============================================================================
# Plot the 3D grid (new version)
# =============================================================================

    def plot(self, scale='log', percent = 98):
        """
        Plot the grid. Do it in a good way, with labels, and with proper orientation.
        
        Optional parameters
        -----
        
        scale:
            The scaling method to use. Can be `log` or `linear`.
            
        """
        
        # For the order, we always want (0,0) to be lower left when we plot.
        
        origin = 'lower'                               # Does imshow() plot start w/ (0,0) at lower-left, or upper-left?
    
        do_stretch_log = ('log' in scale)
        
        # Define the axes. DK says that the array is delivered in order [x, y, z], which are same as Mark's coord frame.
        # That means that if I sum in the '0' direction, I will have a plot in Y and Z.
        
        axes           = ['X',     'Y',     'Z']    
        num_axis       = {'X' : 0, 'Y' : 1, 'Z' : 2}
        
        # Now, make a dictionary to show us what the axes will be of the output image after doing np.sum().
        # The dictionary here means: 
        #     If we sum along the X axis (0), and we plot the result, what will be on vertical axis of the imshow() plot.
        #     In this case, first remaining axis (Y) is on vertical, and second (Z) is on horizontal.
        #     That is how imshow() works.
        
        axes_vertical  = {'X':'Y', 'Y':'X', 'Z':'X'}   
        axes_horizontal= {'X':'Z', 'Y':'Z', 'Z':'Y'}
        axes_transpose = {'X':True,'Y':True,'Z':True}
        view           = {'X':'Side', 'Y':'Front', 'Z':'Top'}  # "If summed along Z axis, this is a Top view", etc.
        
        # Set the font size
        
        hbt.fontsize(10)
        fontsize_axes = 15
        width_colorbar = 10
                
        halfwidth_km = self.km_per_cell_x * hbt.sizex(self.density) / 2
        extent = [-halfwidth_km, halfwidth_km, -halfwidth_km, halfwidth_km]  # Make calibrated labels for X and Y axes
    
        i = 1  # Index of which plot we are on
        
        for axis in axes:

            plt.subplot(1,3,i)

            # Plot the individual image. Start it with origin at lower-left corner.

            img = np.sum(self.density, axis=num_axis[axis])  # Make the flattened image

            if do_stretch_log:
                img_stretch = stretch_hbt(img)
            else:
                img_stretch = astropy.visualization.PercentileInterval(percent)(img)
            
            # Create the colorbar, and superimpose it on the image
            
            colorbar = hbt.frange(np.amin(img_stretch), np.amax(img_stretch), hbt.sizey(img_stretch))
            
            # Manually draw the colorbar onto the plot
            
            if axes_transpose[axis]:
                for j in range(width_colorbar):
                    img_stretch[-j,:] = colorbar  # There is probably a better way to do this?
        
            else:
                for j in range(width_colorbar):
                    img_stretch[:,-j] = colorbar  # There is probably a better way to do this?
                
            # Display the image.
            # The image is displayed in exactly the same orientation as if I print it, with the exception
            # that the origin={lower | upper} keyword can flip it vertically.
            # When accessing the array elements, they are in order img[y, x] -- which is opposite IDL.
            
            if axes_transpose[axis]:
                plt.imshow(np.transpose(img_stretch), extent=extent, origin=origin)
            else:
                plt.imshow(             img_stretch,  extent=extent, origin=origin)
            
            # Create the labels for the colorbar, and individually place them
            
            num_ticks_colorbar = 5 # Number of vertical value to put on our colorbar
            
            for j in range(num_ticks_colorbar):
                val = stretch_hbt_invert(hbt.frange(np.amin(img_stretch), np.amax(img_stretch), num_ticks_colorbar)[j])
                val = round(val)  # Convert to zero if it is very close
                xval = 0.65*np.max(extent)
                yrange = np.max(extent)-np.min(extent)
                yval = np.min(extent) + 0.02*yrange + (j/(num_ticks_colorbar-1) * 0.92*yrange)
                if do_stretch_log:
                    plt.text(xval, yval, f'{val:.1e}', color = 'white') # Label the colorbar, if log stretch only.
            
            # Label the axes and the plot
            
            if axes_transpose[axis]:
                plt.title(f'Summed along {axis} ({view[axis]})', fontsize=fontsize_axes)
                plt.xlabel(axes_vertical[axis] + ' [km]', fontsize=fontsize_axes)
                plt.ylabel(axes_horizontal[axis] + ' [km]', fontsize=fontsize_axes)
                plt.tight_layout()

            else:
                plt.title(f'Summed along {axis}', fontsize=fontsize_axes)
                plt.ylabel(axes_vertical[axis] + ' [km]', fontsize=fontsize_axes)
                plt.xlabel(axes_horizontal[axis] + ' [km]', fontsize=fontsize_axes)
                plt.tight_layout()
            
            i+=1
            
        plt.show()
        
# =============================================================================
# Print some info about the run
# =============================================================================

    def print_info(self):

        print(f'Name:       {self.name}')
        print(f'Beta:       {self.beta}')
        print(f'Duration:   {self.duration}')
        print(f'Speed:      {self.speed}')
        print(f'Timestep:   {self.time_step}')
        print(f'Body ID:    {self.body_id}')
        print(f'Subset:     {self.subset}')
        print(f'# grains:   {self.grains}')
        print(f'State file: {self.state_file}')
        print(f'Obj   file: {self.obj_file}')
        
        return

# =============================================================================
# Write a single file out to an xyz grid, for use in a point cloud mesh
# =============================================================================

    def write_file_xyz(self):
        """
        Write the current grid into a file as an XYZ list of points, rather than a density per grid.
        This might be useful for rendering in other packages, such as MeshLab.
        
        MeshLab by default uses format "X; Y; Z"
        It can also use                "X; Y; Z; R; G; B", where RGB are 0 .. 255 
        
        Many other formats too, but these are easiest.
        
        I would prefer to do rendering within Python, but I could not get IPyvolume, etc. working.
        
        """

        file_out_xyz = self.file_header_full.replace('header.txt', 'xyz.txt')
        
        f = open(file_out_xyz, 'w')
        
        (nx, ny, nz) = np.shape(self.density)
        for i_x in range(nx):
            for i_y in range(ny):
                for i_z in range(nz):
                      j = 0
                      if self.density[i_x, i_y, i_z]:
                          f.write(f'{i_x/nx - 0.5:0.3f}; {i_y/ny - 0.5:0.3f}; {i_z/nz - 0.5:0.3f}; {i_y}; {i_z}; 0\n')
        f.write('1.00')  # I don't know what this is. But the file seems to require it.
        f.close()
        print(f'Wrote: {file_out_xyz}')                 
        
# =============================================================================
# END OF CLASS DEFINITONS
# =============================================================================
        
# =============================================================================
# Make a plot showing all of the runs in the grid. 
# =============================================================================

def plot_flattened_grids(arr):
    
    """
    Plot all of the grids in an array of grids. This is just to give a quick 'family portrait' overview.
    They are plotted in sequential order, but no attempt is made to list by beta, etc.

    parameters
    -----
    
    arr: 
        3D array, with each plane to be plotted.

    """
    
#    stretch_percent = 98  # Since bg is mostly black, this requires a slightly different stretch than normal!
#    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
 
    num = hbt.sizex(arr) # Number of runs in the array
    num_grid_x = math.ceil(math.sqrt(num))
    num_grid_y = num_grid_x
    dx_pix_grid = hbt.sizey(arr)
    dy_pix_grid = dx_pix_grid
    
    img = np.zeros((num_grid_x*dx_pix_grid, num_grid_y*dy_pix_grid))
    i = 0  # Row
    j = 0  # Column
    for k in range(len(runs_full)):
        
        # For each one, copy the grid, normalize it, and put it into out
        
        img[i*dx_pix_grid:(i+1)*dx_pix_grid, j*dy_pix_grid:(j+1)*dy_pix_grid] = \
            hbt.normalize( density_flattened_arr[k,:,:] )                         
        i += 1
        
        if (i > (num_grid_x-1)):
            i = 0
            j += 1
            
    plt.imshow(stretch_hbt(img))
    plt.show()    
    
# =============================================================================
# Make a plot showing all of the runs in the grid. 
# *** This is the same as the function above -- just uses input in table, rather than array! ***
# =============================================================================

def plot_flattened_grids_table(t, stretch_percent=96, file_out = None):
    
    """
    Plot all of the grids in an array of grids. This is just to give a quick 'family portrait' overview.
    They are plotted in sequential order, but no attempt is made to list by beta, etc.
    
    parameters
    -----
    
    t: 
        Astropy table, with column 'img_2d' to be plotted.
    """
 
    stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
    num = hbt.sizex(t['img_2d']) # Number of runs in the array
    num_grid_x = math.ceil(math.sqrt(num))
    num_grid_y = num_grid_x
    dx_pix_grid = hbt.sizey(t['img_2d'][0])
    dy_pix_grid = dx_pix_grid
    
    img = np.zeros((num_grid_x*dx_pix_grid, num_grid_y*dy_pix_grid))
    i = 0  # Row
    j = 0  # Column
    for k in range(num):        

        # For each one, copy the grid, normalize it, and put it into out

        img[j*dy_pix_grid:(j+1)*dy_pix_grid, i*dx_pix_grid:(i+1)*dx_pix_grid] = \
            hbt.normalize( t['img_2d'][k] )                         
        i += 1
        
        if (i > (num_grid_x-1)):  # If we fill up a row, go to the next one.
            i = 0
            j += 1
            
    plt.imshow(stretch(img), origin='lower')
    i = 0
    j = 0
    
    fontsize = 8   # Test font size
    
    for k in range(num):
        s = f"pv={t[k]['albedo']}, q={t[k]['q']}, rho={t[k]['rho']:.2f}, v={t[k]['speed']}"
        plt.text(2 + i*dx_pix_grid, 2 + j*dy_pix_grid, s, color='white', fontsize = fontsize)
        i += 1
        if (i > (num_grid_x-1)):  # If we fill up a row, go to the next one.
            i = 0
            j += 1
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    if file_out:
        plt.savefig(file_out)
        print(f'Wrote: {file_out}')
    plt.show()    

# =============================================================================
# Define a good stretch to use. Same idea as the Astropy Visualization stretch, just better for these data.
# =============================================================================

def stretch_hbt(arr):
    return np.log(10 + arr)

def stretch_hbt_invert(arr):
    return np.exp(arr) - 10
        
# =============================================================================
# Simple test routine if file is run directly
# =============================================================================

if (__name__ == '__main__'):
    
    hbt.figsize(13,13)
    
    plt.set_cmap('plasma')
    
    do_out_xyz = False   # Do we write a file that MeshLab.app can read?
    
    dir_dph2 = '/Users/throop/data/ORT2/hamilton/deliveries/'
    dir_dph3 = '/Users/throop/data/ORT3/hamilton/deliveries/'
    dir_dph4 = '/Users/throop/data/ORT4/hamilton/deliveries/'
    
    dir_dph = dir_dph4
    
    dir_dk  = '/Users/throop/data/ORT_Nov18/kaufmann/deliveries/'
    
    
    # runs_full = ['ort4_bc3_10cbr2_dph_3moon/ort4-0003/y3.0/beta1.0e-06/subset00',
    #              'ort4_bc3_10cbr2_dph_3moon/ort4-0003/y3.0/beta1.0e-06/subset01',
    #              'ort4_bc3_10cbr2_dph_3moon/ort4-0003/y3.0/beta1.0e-06/subset02',
    #              ]
    
    runs_full = ['chr3_sunflower3.5k/chr3-0002/y2.2/beta1.0e+00/subset00',
                 'chr3_sunflower3.5k/chr3-0002/y2.2/beta1.0e+00/subset01',
                 'chr3_sunflower3.5k/chr3-0002/y2.2/beta1.0e+00/subset02',
                 ]

    for run in runs_full:
        
        if ('syntrue' in run) or ('chr3' in run):
            dir = dir_dk
            
        if ('DPH' in run) or ('syn_ort' in run) or ('dph' in run):
            dir = dir_dph

        ring = nh_ort_track3_read(run, dir_base = dir)
        ring.print_info()
        
        # Now make plots of the ring. Make two plots 
        
        plot_lin = False
        plot_log = True
        
        if plot_lin:
            ring.plot(scale='lin')
            
        if plot_log:
            ring.plot(scale='log')
            
        print('\n-----\n')
        
        if do_out_xyz:
            ring.write_file_xyz()
        