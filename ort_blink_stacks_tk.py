#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 09:58:37 2018

This is my image stack blinker for the MU69 ORT. It is an animated Tk tool that allows the user to scan through 
stacked images.

Hit 'h' for help.

Keypresses like this are accepted when cursor is on the GUI. Sometimes, for entering data such as a new timestep, 
text must be entered into the main Python window, not the GUI.

@author: throop
"""

# General python imports

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
from   subprocess import call
import warnings
import pdb
import os.path
import os
import subprocess

import astropy
import astropy.visualization
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
import time
from   importlib import reload            # So I can do reload(module)

import random

import re # Regexp
import pickle # For load/save
import time

# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
import tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

from   astropy.stats import sigma_clip

from   image_stack import image_stack

from get_radial_profile_circular import get_radial_profile_circular

# HBT imports

import hbt

class App:

##########
# INIT CLASS
##########

    def __init__(self, master, size_window):

        self.master = master   # This is the handle to the main Tk widget. I have to use it occasionally to 
                               # set up event handlers, so grab it and save it.

        self.size_window = size_window # Save the size of the whole Tk window, in pixels.
        
        # Open the image stack
        
#        self.stretch_percent = 90    
#        self.stretch = astropy.visualization.PercentileInterval(self.stretch_percent) # PI(90) scales to 5th..95th %ile.
#        
        name_ort = 'ORT4'
#        name_ort = 'ORT2_OPNAV'
        
        if (name_ort == 'ORT1'):
            self.reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
            #        self.reqids_haz  = ['K1LR_HAZ03', 'K1LR_HAZ01', 'K1LR_HAZ02']
            self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'        
            self.dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'

        if (name_ort == 'ORT2'):
            self.reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
#            self.reqids_haz  = ['K1LR_HAZ03', 'K1LR_HAZ01', 'K1LR_HAZ02']
            self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'        
            self.dir_data    = '/Users/throop/Data/ORT2/throop/backplaned/'

        if (name_ort == 'ORT3'):
            self.reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
#            self.reqids_haz  = ['K1LR_HAZ03', 'K1LR_HAZ01', 'K1LR_HAZ02']
            self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'        
            self.dir_data    = '/Users/throop/Data/ORT3/throop/backplaned/'

        if (name_ort == 'ORT4'):
            self.reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03']
            self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'        
            self.dir_data    = '/Users/throop/Data/ORT4/throop/backplaned/'
            
        if (name_ort == 'ORT2_OPNAV'):
            self.dir_data    = '/Users/throop/Data/ORT2/throop/backplaned/'
            dirs = glob.glob(self.dir_data + '/*LR_OPNAV*')         # Manually construct a list of all the OPNAV dirs
            self.reqids_haz = []
            for dir_i in dirs:
                self.reqids_haz.append(os.path.basename(dir_i))
            self.reqids_haz = sorted(self.reqids_haz)    
            self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'    
            
        # Set the edge padding large enough s.t. all output stacks will be the same size.
        # This value is easy to compute: loop over all stacks, and take max of stack.calc_padding()[0]
        
        self.padding     = 61 # Amount to pad the images by. This is the same as the max drift btwn all images in stacks
        self.zoom        = 2  # Sub-pixel zoom to apply when shifting images. 1 for testing; 4 for production.
        self.num_image   = 0  # Which stack number to start on.
        self.zoom_screen = 1  # 'Screen zoom' amount to apply. This can be changed interactively.
        
        self.is_blink    = False  # Blinking mode is turned off by default
        self.dt_blink    = 300    # Blink time in ms
        
        # Start up SPICE if needed
        
        if (sp.ktotal('ALL') == 0):
            sp.furnsh('kernels_kem_prime.tm')
            
        # Set the RA/Dec of MU69. We could look this up from SPICE but it changes slowly, so just keep it fixed for now.
        
        self.radec_mu69 = (4.794979838984583, -0.3641418801015417)
        
        # Set the CA time. Roughly doing this is fine.
        
        self.et_ca = sp.utc2et('2019 1 Jan 05:33:00')
        
        # Boolean. For the current image, do we subtract the field frame, or not?
        
        self.do_subtract = True

        hbt.set_fontsize(20)

        # Set the stretch range, for imshow. These values are mapped to black and white, respectively.
        
        self.vmin_diff = -1   # Range for subtracted images
        self.vmax_diff =  2
        
        self.vmin_raw = -1    # Range for raw images (non-subtracted)
        self.vmax_raw = 1000
        
# Restore the stacks directly from archived pickle file, if it exists
        
        self.file_save = os.path.join(self.dir_data, 
                                      f'stacks_blink_{name_ort}_n{len(self.reqids_haz)}_z{self.zoom}.pkl')
        
        if os.path.isfile(self.file_save):
            self.restore()
        else:

# If no pickle file, load the stacks from raw images and re-align them
            
            # Load and stack the field images
    
            print("Stacking field images")        
            self.stack_field = image_stack(os.path.join(self.dir_data, self.reqid_field))    # The individual stack
            self.stack_field.align(method = 'wcs', center = self.radec_mu69)
            (self.img_field, self.wcs_field)  =\
                self.stack_field.flatten(zoom=self.zoom, padding=self.padding) # Save the stacked image and WCS
        
            # Load and stack the Hazard images
            
            self.img_haz   = {} # Output dictionary for the stacked images
            self.stack_haz = {} # Output dictionary for the stacks themselves
            self.wcs_haz   = {} # Output dictionary for WCS for the stacks
            
            for reqid in self.reqids_haz:
                self.stack_haz[reqid] = image_stack(os.path.join(self.dir_data, reqid))    # The individual stack
                self.stack_haz[reqid].align(method = 'wcs', center = self.radec_mu69)
                (self.img_haz[reqid], self.wcs_haz[reqid])  =\
                    self.stack_haz[reqid].flatten(zoom=self.zoom, padding=self.padding) 
                # Put them in a dictionary

            # Save the stacks to a pickle file, if requested
            
            yn = input("Save stacks to a pickle file? ")
            if ('y' in yn):
                self.save()
                
# Set the sizes of the plots -- e.g., (15,15) = large square
        
        figsize_image = (15,15)
        
        self.fig1 = Figure(figsize = figsize_image)    # <- this is in dx, dy... which is opposite from array order!

        self.ax1 = self.fig1.add_subplot(1,1,1, 
                                    label = 'Image') # Return the axes
        plt.set_cmap('Greys_r')
        
        self.canvas1 = FigureCanvasTkAgg(self.fig1,master=master)
        self.canvas1.show()
        
# Put objects into appropriate grid positions

        self.canvas1.get_tk_widget().grid(row=1, column=1, rowspan = 1)
        
# Define some keyboard shortcuts for the GUI
# These functions must be defined as event handlers, meaning they take two arguments (self and event), not just one.

        master.bind('q',       self.quit_e)
        master.bind('<space>', self.toggle_subtract_e)
        master.bind('=',       self.prev_e)
        master.bind('-',       self.next_e)
        master.bind('h',       self.help_e)
        master.bind('?',       self.help_e)
        master.bind('<Left>',  self.prev_e)
        master.bind('<Right>', self.next_e)
        master.bind('s',       self.stretch_e)
        master.bind('b',       self.blink_e)
        master.bind('t',       self.blink_set_time_e)
        master.bind('#',       self.blink_set_sequence_e)
        master.bind('z',       self.zoom_screen_up_e)
        master.bind('Z',       self.zoom_screen_down_e)
        master.bind('x',       self.clear_current_objects_e)
        master.bind('X',       self.clear_all_objects_e)
        
        master.bind('=',       self.scale_max_up_e)
        master.bind('+',       self.scale_max_down_e)
        master.bind('-',       self.scale_min_up_e)
        master.bind('_',       self.scale_min_down_e)
        
        master.bind('S',       self.save_output_e)
        
        self.canvas1.get_tk_widget().bind("<Button 1>", self.click_e)        
        
# Set the initial image index
        
        self.reqid_haz = self.reqids_haz[self.num_image]  # Set it to 'K1LR_HAZ00', for instance.

# Initialize the list of found objects for each stack
# There is a list of objects for each individual stack (ie, for each frame in the blink)

        self.list_objects = {}
        
        for reqid_i in self.reqids_haz:
            self.list_objects[reqid_i] = []  # Each entry here will be something like [(x, y, dn), (x, y, dn)] 

# Initialize a set of matplotlib 'line' objects for the image.
# These correspond to the 'objects' above, which are really just points            
            
        self.list_lines = {}

        for reqid_i in self.reqids_haz:
            self.list_lines[reqid_i] = []  # Each entry here will be a list of plot objects, of type 'line' 
                    
# Set a list of frame numbers to animate. For default, do them all.

        self.list_index_blink = hbt.frange(0, len(self.reqids_haz)-1) # List of indices ( [1, 2, 3] )
        self.list_index_blink_str = ' '.join(np.array(self.list_index_blink).astype(str)) # Make into string ('1 2 3')
        self.index_blink = 0      # where in the list of indices do we start? Current index.     
        
# Plot the image
        
        self.plot()


# =============================================================================
# Restore all stacks from a single saved pickle array
# =============================================================================

    def restore(self):
        
        """
        Load state from a pickle file. No arguments.
        """

        print("Reading: " + self.file_save)           
        lun = open(self.file_save, 'rb')
        (self.stack_field, self.img_field, self.stack_haz, self.img_haz, self.wcs_haz) = pickle.load(lun)
        lun.close()

        return
    
# =============================================================================
# Save current state into a pickle file.
# =============================================================================

    def save(self):
        
        """
        Save stacks into a .pkl file. No arguments.
        """

        lun = open(self.file_save, 'wb')
        pickle.dump((self.stack_field, self.img_field, self.stack_haz, self.img_haz, self.wcs_haz), lun)
        lun.close()
        print("Wrote: " + self.file_save)
        
        return
        
# =============================================================================
# Plot the current image, updating axes etc.
# =============================================================================
        
    def plot(self):

#        print("in self.plot()")
        
        # Load the current image
                
        img_haz = self.img_haz[self.reqids_haz[self.num_image]]
        
        # Calculate and apply the 'screen zoom'
        
        dx = hbt.sizex(img_haz)
        dx_zoomed = dx / self.zoom_screen
        min = int(dx/2 - dx_zoomed/2)
        max = int(dx/2 + dx_zoomed/2)
        
        origin = 'upper'  # Set (0,0) to be uper-left corner, which is default. However, 
        
        # Clear the existing image. Otherwise, we end up over-plotting, and the speed drops every time plot() called.
        
        self.ax1.clear()
        
        # Subtract the bg image, if requested, and display it
        
        # Note that we don't use xlim / ylim here. We just crop the array. I'm not sure about that -- I think 
        # maybe we should use xlim / ylim, so as to keep the XY position of MU69 constant. It's just going to 
        # make the math a lot easier.

        if (self.do_subtract):
            im_disp = img_haz - self.img_field
            self.plt1 = self.ax1.imshow(im_disp, interpolation=None, vmin=self.vmin_diff, vmax=self.vmax_diff,
                                        origin = origin)
            self.ax1.set_ylim((min, max))
            self.ax1.set_xlim((min, max))
            
        else:
            self.plt1 = self.ax1.imshow(img_haz, interpolation=None, 
                                        vmin=self.vmin_raw, vmax=self.vmax_raw,
                                        origin = origin)
            
        # Set the axis limits
        
        self.ax1.set_ylim((min, max))
        self.ax1.set_xlim((min, max))

#        print(f'Min={min}, max={max}')
        
        # Plot any astrometric points that may exist
        
        self.plot_objects()
        
        # Set the title, etc.
        
        self.ax1.set_title(self.make_title())
            
        # Make it as compact as possible
        
        self.fig1.tight_layout()
        
        self.canvas1.show()


    def make_title(self):
        """
        Return the title for a plot
        """
        
#        str_et = 
        return 'K{:3.1f}d, {}, {}/{}, zoom = {},{}, vscale=[{:4.2f},{:4.2f}]'.format(
                       (self.stack_haz[self.reqid_haz].et - self.et_ca)/86400,
                       self.reqid_haz,
                       self.num_image, len(self.reqids_haz), self.zoom, self.zoom_screen,
                       self.vmin_diff, self.vmax_diff)
    
# =============================================================================
# Update just the data in the current image, leaving all the rest untouched.
# =============================================================================
        
    def replot(self):

#        print("in self.replot()")

        # Remove all of the current points
        
        self.unplot_objects()
        
        # Load the current image

        # Figure out key for newest image to plot
        self.reqid_haz = self.reqids_haz[self.num_image]  # Set it to 'K1LR_HAZ00', for instance.

#        print(f'Plotting num_image={self.num_image} : self.reqid_haz={self.reqid_haz}')            
        
        # Load the image
        img = self.img_haz[self.reqid_haz]
        
#        self.img_haz = self.stack_haz.image_single(self.num_image, padding = self.padding, zoom = self.zoom)

#        dx = hbt.sizex(img)
#        dx_zoomed = dx / self.zoom_screen
#        min = int(dx/2 - dx_zoomed/2)
#        max = int(dx/2 + dx_zoomed/2)
        
        # Subtract the bg image, if requested, and display it
        
        if (self.do_subtract):
#            im_disp = img[ - self.img_field[min:max, min:max]
            im_disp = img - self.img_field
            self.plt1.set_data(im_disp)
        else:
#            self.plt1.set_data(img[min:max, min:max]) # Like imshow, but faster
            self.plt1.set_data(img) # Like imshow, but faster

        # Set the title. This gets updated when using draw() or show() just fine.
                
        self.ax1.set_title(self.make_title())
        

        
        # Draw the photometric points for this frame. 
        
        self.plot_objects()
        
        self.canvas1.draw()
        self.canvas1.show()  # Q: Do I need this:? A: Yes, even when using set_data().
        
# =============================================================================
# Key: Help
# =============================================================================

    def help_e(self, event):
        print('------------------------------')
        print('Help!')
        print('<space>      : Toggle field subtraction on/off')
        print('Z / z        : Zoom in/out' )
        print('<left> or =  : Previous')
        print('<right> or - : Next')
        print('-            : Next')
        print('h or ?       : Help')
        print('q            : Quit')
        print('s            : Change stretch')
        print('b            : Toggle blinking on/off')
        print('t            : Change blink time')
        print('#            ; Change blink sequence')
        print('+-           : Change vertical scaling')
        print('+- shifted   : Change vertical scaling')
        print('x            : Clear all points from current frame')
        print('X            : Clear all points from all     frames')
        print('S            : Save astrometric output to a file')
        
        print('------------------------------')
        print('CURRENT STATUS:')
        print(f'Blink sequence: {self.list_index_blink_str}')
        print(f'Blink status: {self.is_blink}, dt = {self.dt_blink}')
        print(f'Zoom:        {self.zoom}')
        print(f'Zoom_screen: {self.zoom_screen}')
        print(f'v1, v2 diff: {self.vmin_diff},  {self.vmax_diff}')
        print(f'v1, v2 raw:  {self.vmin_raw},   {self.vmax_raw}')
        
# =============================================================================
# Key: Next image
# =============================================================================

    def next_e(self, event):
        self.num_image_prior = self.num_image  # Save it so we can blink to it.
        self.num_image += 1
        self.num_image = np.clip(self.num_image, 0, len(self.reqids_haz)-1)
        print(f'Next: num_image = {self.num_image}')
        print(f'Next: stack = {self.reqids_haz[self.num_image]}')
        self.replot()

# =============================================================================
# Key: Previous image
# =============================================================================
        
    def prev_e(self, event):
        """
        Go to previous image. This is *not* the previous image in blink sequence.
        """
        
        self.num_image_prior = self.num_image
        self.num_image -=1
        self.num_image = np.clip(self.num_image, 0, len(self.reqids_haz)-1)
        
        print(f'Prev: num_image = {self.num_image}')
        print(f'Prev: stack = {self.reqids_haz[self.num_image]}')
        self.replot()


# =============================================================================
# Key: Clear current objects
# =============================================================================

    def clear_all_objects_e(self, event):

        s = input(f'Really clear all objects from all frames? ')

        if ('y' in s):
            for reqid_i in self.reqids_haz:
                self.list_objects[reqid_i] = []
            
            self.unplot_objects()
        
# =============================================================================
# Key: Clear current objects
# =============================================================================

    def clear_current_objects_e(self, event):

        s = input(f'Really clear all objects from this frame? ')

        if ('y' in s):    
            self.list_objects[self.reqid_haz] = []            
            self.unplot_objects()        
        
        self.plot()
        
# =============================================================================
# Key: Change Stretch
# =============================================================================
        
    def stretch_e(self, event):
        """
        Change the min and max stretch values. We scale linearly between these.
        These are set separately for the case of subtracting the bg, and not.
        """

        if self.do_subtract:
            s = input(f'Enter new diff stretch range ({self.vmin_diff} {self.vmax_diff}): ')
            if s:
                self.vmin_diff = int(s.strip().split(' ')[0])  # Mapped to black
                self.vmax_diff = int(s.strip().split(' ')[1])  # Mapped to white 
        else:
            s = input(f'Enter new raw stretch range ({self.vmin_raw} {self.vmax_raw}): ')
            if s:
                self.vmin_raw = int(s.strip().split(' ')[0])  # Mapped to black
                self.vmax_raw = int(s.strip().split(' ')[1])  # Mapped to white 
        
        self.plot()   # We have to do plot and not replot, since the latter doesn't take vmin/vmax args.

# =============================================================================
# Key: zoom screen Up
# =============================================================================

    def zoom_screen_up_e(self, event):
        self.zoom_screen += 1
        print(self.zoom_screen)
        self.plot()

# =============================================================================
# Key: zoom screen Down
# =============================================================================

    def zoom_screen_down_e(self, event):
        self.zoom_screen -= 1
        if (self.zoom_screen < 1):
            self.zoom_screen = 1
        print(self.zoom_screen)
        self.plot()

# =============================================================================
# Mouse click handler
# =============================================================================

    def click_e(self, event):

        (x, y) = (event.x, event.y)                             # Event coords. These are pixels from upperleft.
        
        y = self.size_window[1]-y                               # Event.xy returns coords from upper left.
                                                                # But transData assumes display coords from lower left.
                                                                # So, need to invert the y axis. Not sure why, but
                                                                # this correction does it properly.
        
        tr_data_to_display   = self.ax1.transData               # Data = values of xlim, ylim
        tr_display_to_data   = self.ax1.transData.inverted()
        
        tr_axes_to_display = self.ax1.transAxes                 # Axes = [0,1] within the plot data window itself.
        tr_display_to_axes = self.ax1.transAxes.inverted()

        tr_fig_to_display  = self.fig1.transFigure              # Figure = [0,1] in the figure, incl borders, etc.
        tr_display_to_fig  = self.fig1.transFigure.inverted()
        
        (x_data, y_data) = tr_display_to_data.transform((x,y))  # *** This provides x+y values in image pixels. ***
        (x_axes, y_axes) = tr_display_to_axes.transform((x,y))
        (x_fig, y_fig) = tr_display_to_fig.transform((x,y))
        
        
#        print(f'Raw: X={x}, Y={y}')
        print(f'Display-to-data: X={x_data}, Y={y_data}') # This works 100% right for x!
                                                          # For y, it is inverted, and has an offset, w/ scale correct. 

#        print(f'Display-to-axes: X={x_axes}, Y={y_axes}') # This is 100% correct. 
#                                                          # (0,0) = upper-left corner of image itself.
#                                                          # (1,1) = lower-right corner of image.
#        
#        print(f'Display-to-fig: X={x_fig}, Y={y_fig}')

        wcs = self.wcs_haz[self.reqid_haz]  # Get WCS of the current stack being displayed. Others should be the same.
        radec = wcs.wcs_pix2world(x_data, y_data, 0)
        print(f'RA = {radec[0]}, Dec = {radec[1]} / x = {x_data}, Y = {y_data}')
#        self.ax1.plot([x_data], [y_data], marker = 'o', color = 'red', ms=5)
#        self.canvas1.show()
                
        print(f'State for this event: {event.state}')   # Check if it was a shift-click, double-click, etc.
                                                        # Shift = 1. Cmd = 8. Alt = 16. Single = 0.
        
        # Now add this point to the list of points for this stack
        
        self.list_objects[self.reqid_haz].append((x_data, y_data, 0))
        
        self.plot()  # Replot the entire window, which will call plot_objects()

# =============================================================================
# Key: Blink On/Off
# =============================================================================

    def blink_e(self, event):
        
        if (self.is_blink):
            self.is_blink = False
            print("Blinking now off")
        
        else:
            self.is_blink = True
            print("Blinking now on")
            self.show_next_frame()            # To turn on animation, set the flag, and call func to display next frame

# =============================================================================
# Plot all of the astrometric points on this stack
# =============================================================================

    def plot_objects(self):
        
        self.list_lines[self.reqid_haz] = []
        
        for i,object in enumerate(self.list_objects[self.reqid_haz]):
#            print(f'plotting object {i}')
#            print(f'object = {object}')
            line = self.ax1.plot(object[0], object[1], marker = 'o', color = 'red', markersize=4)[0]
            
            # Add this 'line' object to the list of lines in this plane
            
            self.list_lines[self.reqid_haz].append(line)

# =============================================================================
# Remove the plotted points from the current image.
# =============================================================================

    def unplot_objects(self):

        # As needed, unplot these objects
        
        lines = self.list_lines[self.reqid_haz]
        
        for i,line in enumerate(lines):
#            print(f'Removing line {i} = {line}')
            line.remove()
        
        return
        
# =============================================================================
# Show next frame -- animation
# =============================================================================

    def show_next_frame(self):
        
        """
        This function when called advances to the next frame in the animation.
        It then sets an 'idle event handler' with a timeout, so that it gets called
        again automatically by Tk. Other than setting this event handler, the animation
        is entirely done in the background, and the app responds fully to all events.
        """
        
        self.index_blink += 1  # Go to the next item in the sequence list
        if self.index_blink == len(self.list_index_blink):
            self.index_blink = 0  # If we're at the end, go back to start
            
        self.num_image = self.list_index_blink[self.index_blink]
        if (self.num_image >= len(self.reqids_haz)):
            self.num_image = 0
#        print(f"Animating frame {self.num_image}")    
        self.replot()    
        if (self.is_blink):
            self.master.after(self.dt_blink, self.show_next_frame)
        
# =============================================================================
# Key: Toggle bg subtraction
# =============================================================================
        
    def toggle_subtract_e(self, event):
        self.do_subtract = not(self.do_subtract)
        self.plot()        
        
# =============================================================================
# Key: Change blink time
# =============================================================================
        
    def blink_set_time_e(self, event):

        dt = input(f'Enter blink time in ms ({self.dt_blink}): ')
        if (dt):                       # If user didn't just hit <cr>
            self.dt_blink = int(dt)

        print(f'dt = {dt} ms')
        
# =============================================================================
# Key: Change blink sequence
# =============================================================================
        
    def blink_set_sequence_e(self, event):

        s = input(f'Enter frame list to animate ({self.list_index_blink_str}): ')
        if (s):
            self.list_index_blink = list(np.array(s.strip().split(' ')).astype(int))
                      
        self.list_index_blink_str = ' '.join(np.array(self.list_index_blink).astype(str))

        self.index_blink = 0           # Start at the initial animation frame
        
        print(f'Blinking frames {self.list_index_blink}')
        
# =============================================================================
# Keys: Change scaling (four different choices)
# =============================================================================
        
    def scale_max_up_e(self, event):
        if (self.do_subtract):
            self.vmax_diff *= 1.1
        else:
            self.vmax_raw *= 1.1            
        self.plot()

    def scale_max_down_e(self, event):
        if (self.do_subtract):
            self.vmax_diff /= 1.1
        else:
            self.vmax_raw /= 1.1            
        self.plot()
        
    def scale_min_up_e(self, event):
        if (self.do_subtract):
            self.vmin_diff *= 1.1
        else:
            self.vmin_raw *= 1.1            
        self.plot()

    def scale_min_down_e(self, event):
        if (self.do_subtract):
            self.vmin_diff /= 1.1
        else:
            self.vmin_raw /= 1.1            
        self.plot()


# =============================================================================
# Key: Save Output to a file
# =============================================================================
        
    def save_output_e(self, event):
                
        """
        Write output to a file, using API as described in
        https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Detection-to-Object
        """
        
        id_user = 7  # HBT personal ID, as e
        inits_user = 'hbt'
        
        for reqid_i in self.reqids_haz:
            for i,pt_i in enumerate(self.list_objects[reqid_i]):
                name_object = f'{id_user}{i:03}'
                wcs = self.wcs_haz[reqid_i]  # Get the proper wcs
                radec = wcs.pix2world(self.list_objects[reqid_i][i])
                ra = radec[0]
                dec = radec[1]
                
# =============================================================================
# Key: Quit
# =============================================================================
        
    def quit_e(self, event):
        self.quit()        

# =============================================================================
# Quit app
# =============================================================================
    
    def quit(self):
        root.destroy() # Q: How does root. get automatically imported here? Not sure.
        root.quit()        
                
###########
# Now start the main app
###########

size_window = (1080, 1080) 

# Start up the widget

root = tkinter.Tk()
app  = App(root, size_window)

# set the dimensions and position of the window

root.geometry('%dx%d+%d+%d' % (size_window[0], size_window[1], 2, 2)) # I chose this size 
                                                    # experimentally so it would match the figsize().
root.configure(background='#ECECEC')                # ECECEC is the 'default' background
               
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

while True:
    try:
        root.mainloop()
        break
    except UnicodeDecodeError:  # Fix a TK + py3 bug, which causes Unicode errors on events. 
                # http://stackoverflow.com/questions/16995969/inertial-scrolling-in-mac-os-x-with-tkinter-and-python
        pass
    
    
#def other:

# https://stackoverflow.com/questions/292095/polling-the-keyboard-detect-a-keypress-in-python   
#        https://stackoverflow.com/questions/29158220/tkinter-understanding-mainloop
        
    
#def other
#
#    dx = 10
#    arr = hbt.dist_center(dx)
#    num_pts = 5
#    
#    pts = []
#    
#    ax = plt.imshow(arr)
#    for i in range(num_pts):
#        pt = plt.plot([random.random()*dx], [random.random()*dx], marker = 'o', ms = 6, ls = 'none', color = 'red')[0]
#        pts.append(pt)
#
#    for i in range(num_pts):
#        if (random.random() > 0.5):
#            print (f'Removing point {i}')
#            pt = pts[i]
#            pt.remove()
#            
#    plt.show()
    