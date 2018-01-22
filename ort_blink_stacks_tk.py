#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 09:58:37 2018

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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
import time
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save


# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
import tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

from   astropy.stats import sigma_clip

from   image_stack import image_stack

# HBT imports

import hbt

class App:


##########
# INIT CLASS
##########

    def __init__(self, master):

        # Open the image stack
        
        stretch_percent = 90    
        self.stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
        
        self.reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
        self.reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        
        self.dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'
            
        # Set the edge padding large enough s.t. all output stacks will be the same size.
        # This value is easy to compute: loop over all stacks, and take max of stack.calc_padding()[0]
        
        self.padding     = 61 # Amount to pad the images by. This is the same as the max drift btwn all images in stacks
        self.zoom        = 1  # Sub-pixel zoom to apply when shifting images
        self.num_image   = 5  # Which image number to start on.
        self.zoom_screen = 1  # 'Screen zoom' amount to apply. This can be changed interactively.
                
        # Start up SPICE if needed
        
        if (sp.ktotal('ALL') == 0):
            sp.furnsh('kernels_kem.tm')
            
        # Set the RA/Dec of MU69. We could look this up from SPICE but it changes slowly, so just keep it fixed for now.
        
        self.radec_mu69 = (4.794979838984583, -0.3641418801015417)
        
        # Boolean. For the current image, do we subtract the field frame, or not?
        
        self.do_subtract = True
        
        # Load and stack the field images
        
        stack_field = image_stack(os.path.join(self.dir_data, self.reqid_field))
        stack_field.align(method = 'wcs', center = self.radec_mu69)
        self.img_field  = stack_field.flatten(zoom=self.zoom, padding=self.padding)
    
        hbt.figsize((12,12))
        hbt.figsize((5,5))
        hbt.set_fontsize(20)
        
        self.reqid_haz = self.reqids_haz[0]
        
        self.stack_haz = image_stack(os.path.join(self.dir_data, self.reqid_haz))
        self.stack_haz.align(method = 'wcs', center = self.radec_mu69)
        

#        img_haz = stack_haz.image_single(self.num_image, padding = self.padding, zoom = self.zoom)
        
#        for reqid in reqids_haz:
#            stack_haz = image_stack(os.path.join(dir_data, reqid))
#            stack_haz.align(method = 'wcs', center = radec_mu69)
#            img_haz  = stack_haz.flatten(zoom=zoom, padding=padding)

# Set the sizes of the plots -- e.g., (7,7) = large square
        
        figsize_image = (20,20)
        
        self.fig1 = Figure(figsize = figsize_image)    # <- this is in dx, dy... which is opposite from array order!

        self.ax1 = self.fig1.add_subplot(1,1,1, 
                                    xlabel = 'X', ylabel = 'Y', 
                                    label = 'Image') # Return the axes
        plt.set_cmap('Greys_r')
        
        self.canvas1 = FigureCanvasTkAgg(self.fig1,master=master)
        self.canvas1.show()
        
# Set up Plot 2 : Radial / Azimuthal profiles
        
        self.bgcolor = 'red'

        self.fig2 = Figure(figsize = (7,3))
        _ = self.fig2.set_facecolor(self.bgcolor)

        self.ax2 = self.fig2.add_subplot(1,1,1, 
                                    xlabel = 'Radius or Azimuth', ylabel = 'Intensity', 
                                    label = 'Plot') # Return the axes
        
        self.canvas2 = FigureCanvasTkAgg(self.fig2,master=master)
        self.canvas2.show()  

# Put objects into appropriate grid positions

# Column 1
        
        self.canvas1.get_tk_widget().grid(row=1, column=1, rowspan = 1)
        self.canvas2.get_tk_widget().grid(row=2, column=1, rowspan = 1)
        
# Define some keyboard shortcuts for the GUI
# These functions must be defined as event handlers, meaning they take two arguments (self and event), not just one.

        master.bind('q', self.quit_e)
        master.bind('<space>', self.toggle_subtract_e)
        master.bind('=', self.prev_e)
        master.bind('-', self.next_e)
        master.bind('<Left>',  self.prev_e)
        master.bind('<Right>', self.next_e)
        
        master.bind('z', self.zoom_screen_up_e)
        master.bind('Z', self.zoom_screen_down_e)
        
# Plot the image
        
        self.plot()

# =============================================================================
# Plot the current image, updating axes etc.
# =============================================================================
        
    def plot(self):

        # Load the current image
        
        self.img_haz = self.stack_haz.image_single(self.num_image, padding = self.padding, zoom = self.zoom)
        
        # Calculate and apply the 'screen zoom'
        
        dx = hbt.sizex(self.img_haz)
        dx_zoomed = dx / self.zoom_screen
        min = int(dx/2 - dx_zoomed/2)
        max = int(dx/2 + dx_zoomed/2)
        
        # Subtract the bg image, if requested, and display it
        
        if (self.do_subtract):
            im_disp = self.img_haz[min:max, min:max] - self.img_field[min:max, min:max]
            self.plt1 = self.ax1.imshow(self.stretch(im_disp), interpolation=None)
        else:
            self.plt1 = self.ax1.imshow(self.stretch(self.img_haz[min:max, min:max]), interpolation=None)
            
        # Set the title, etc.
        
        self.ax1.set_title('{}, {}/{}, zoom = {}'.format(self.reqid_haz,
                       self.num_image, self.stack_haz.size[0], self.zoom))
#        self.ax1.set_xlim(range)
#        self.ax1.set_ylim(range)
        self.canvas1.show()

#        self.ax2.plot(self.img_haz[500,:])
        self.canvas2.show()
 
# =============================================================================
# Update just the data in the current image, leaving all the rest untouched.
# =============================================================================
        
    def replot(self):

        # Load the current image
        
        self.img_haz = self.stack_haz.image_single(self.num_image, padding = self.padding, zoom = self.zoom)

        dx = hbt.sizex(self.img_haz)
        dx_zoomed = dx / self.zoom_screen
        min = int(dx/2 - dx_zoomed/2)
        max = int(dx/2 + dx_zoomed/2)
        
        # Subtract the bg image, if requested, and display it
        
        if (self.do_subtract):
            im_disp = self.img_haz[min:max, min:max] - self.img_field[min:max, min:max]
            self.plt1.set_data(self.stretch(im_disp))
        else:
            self.plt1.set_data(stretch(im_disp))
                    
        self.canvas1.show()  # Q: Do I need this:? A: Yes, even when using set_data().
        
# =============================================================================
# Key: Next image
# =============================================================================

    def next_e(self, event):
        self.num_image += 1
        self.replot()

# =============================================================================
# Key: Previous image
# =============================================================================
        
    def prev_e(self, event):
        self.num_image -= 1
        self.replot()

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
# Key: Toggle bg subtraction
# =============================================================================
        
    def toggle_subtract_e(self, event):
        self.do_subtract = not(self.do_subtract)
        self.plot()
        
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

# Start up the widget

root = tkinter.Tk()
app  = App(root)

# set the dimensions and position of the window

root.geometry('%dx%d+%d+%d' % (810, 800, 2, 2))
root.configure(background='#ECECEC')                # ECECEC is the 'default' background
               
os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

while True:
    try:
        root.mainloop()
        break
    except UnicodeDecodeError:  # Fix a TK + py3 bug, which causes Unicode errors on events. 
                # http://stackoverflow.com/questions/16995969/inertial-scrolling-in-mac-os-x-with-tkinter-and-python
        pass