# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

@author: throop

"""

# -*- coding: utf-8 -*-
"""
# Program display a GUI of NH J-ring data to navigate and extract it.

# Possible features to be added:
#  o Output images as PNG
#  o Select multiple images and display them as an arrayf
#  o Slider to scale the image brightness / contrast
#  o Subtract median? Subract neighbor image?

# Widget program created 18-May-2016 based on previous text version (same name), and salt_interact.py .

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
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
from tkinter import ttk
from tkinter import messagebox
tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

# HBT imports

import hbt

# Local NH rings imports

from  nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist
#from  nh_jring_load_objectlist      import nh_jring_load_objectlist

from nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image
#from nh_jring_extract_profiles_from_unwrapped  import nh_jring_extract_profiles_from_unwrapped   
from nh_jring_extract_profile_from_unwrapped   import nh_jring_extract_profile_from_unwrapped   

#pdb.set_trace()

# First we define any general-purpose functions, which are not part of the class/module.
# We can move these to a different file at some point.

        
class App:

# Now define the functions that are part of this module (or class, not sure). They are accessed with 
# self.function(), and they have access to all of the variables within self. .	

##########
# INIT CLASS
##########

    def __init__(self, master):

# Set some default values

        self.filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters 

        self.dir_out    = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
        
        dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

        file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel

# Set the sizes of the plots -- e.g., (7,7) = large square

        self.figsize_1 = (7,7) # Setting this larger and larger still leaves whitespace...
        self.figsize_2 = (7,3)
        self.figsize_3 = (7,4)

# Set various default values
        
        option_bg_default   = 'String' # Default backgroiund type. Need to set this longer too.
        entry_bg_default    = '0-10'   # Default polynomial order XXX need to set a longer string length here!
        index_group_default = 7        # Default group to start with
        index_image_default = 32       # Default image number within the group

        self.do_autoextract     = 1             # Flag to extract radial profile when moving to new image. 
                                                # Flag is 1/0, not True/False, as per ttk.
                                                
# Do some general code initialization

        hbt.set_plot_defaults()
        
        self.bgcolor = '#ECECEC'    # ECECEC is the 'default' background for a lot of the ttk widgets.
        
        (r2d, d2r) = (hbt.r2d, hbt.d2r)

        self.stretch_percent = 90    
        self.stretch = astropy.visualization.PercentileInterval(self.stretch_percent) # PI(90) scales to 5th..95th %ile.

# Set some physical parameters

        self.a_ring_inner_km = 127500  # For clarify, we want to plot inward of the real edge
        self.a_ring_outer_km = 129000  # 129000 is slightly inside outer edge -- easy to see
        
# Start up SPICE

        sp.furnsh(file_tm) # Commented out for testing while SPICE was recopying kernel pool.

# Check if there is a pickle save file found. If it is there, go ahead and read it.

        if os.path.exists(self.dir_out + self.filename_save):
            print("Loading file: " + self.dir_out + self.filename_save)
            self.load(verbose=False) # This will load self.t
            t = self.t

        else:    

# Find and process all of the FITS header info for all of the image files
# 't' is our main data table, and has all the info for all the files.
# We want to locate only the 'opnav' files, which are the ones which have correctly navigated WCS info.
        
            t = hbt.get_fits_info_from_files_lorri(dir_images, pattern='opnav')

# Define new columns for the main table. 
# For several of these, we want to stick an array into the table cell. We do this by defining it as a string.
         
            t['dx_opnav']     = 0 # The measured offset between pos_star_cat and pos_star_image
            t['dy_opnav']     = 0
            t['dx_offset']    = 0 # Additional user-requested offset between the two
            t['dy_offset']    = 0
            t['bg_method']    = option_bg_default # None, Previous, Next, Polynomial, Median, "Grp Num Frac Pow", etc.
            t['bg_argument']  = entry_bg_default # A number (if 'Polynomial'). A range (if 'Median')
            t['Comment']      = 'Empty comment'  # Blank -- I'm not sure how to init its length if needed
            t['is_navigated'] = False  # Flag

# Extend the size of these strings. This is really stupid, but Tables pre-allocates width of these fields, and 
# will quietly truncate any longer string we stick in there. So, need to explicitly make them wider.

            t['bg_argument'] = np.array(t['bg_argument'], dtype = 'U100')
            t['bg_method']   = np.array(t['bg_method'],   dtype = 'U100')
            t['Comment']     = np.array(t['Comment'],     dtype = 'U100')
                      
            self.t = t
            
            # Since there was no pickle file read, go ahead and write it out. Some of the other routines
            # rely on having this file available.

            lun = open(self.dir_out + self.filename_save, 'wb')
            t = self.t
            pickle.dump(t, lun)
            lun.close()
            print("Wrote: " + self.dir_out + self.filename_save)
            
        # Get a list of unique observation descriptions (i.e., 'groups')
        
        self.groups = astropy.table.unique(t, keys=(['Desc']))['Desc']
                
        # Set initial location to first file in first group
        
        self.index_group = index_group_default
        self.index_image  = index_image_default
        
        # Extract the one group we are looking at 
        # NB: This creates an entirely new table -- *not* a view into the original table. 
        # If we modify values here, we need to explicitly write the changes back to the original.
        
        self.groupmask = t['Desc'] == self.groups[self.index_group]
        self.t_group = t[self.groupmask]  # 
        
        self.num_images_group = np.size(self.t_group)

# Now set up a few variables for the current image only. These are not stored in the main
# table, because they are regenerated / reloaded every time, and they are potentially very large.

        self.planes             = 0         # Backplanes will be loaded into this variable
        self.has_backplane      = False     # Need to update this every time we load a new image file.
        self.objectlist         = 0
        self.has_objectlist     = False
        self.profile_radius     = np.zeros((1))
        self.profile_azimuth    = np.zeros((1))
        self.bins_radius        = np.zeros((1))
        self.bins_azimuth       = np.zeros((1))
        
        self.image_raw          = np.zeros((1)) # The current image, with *no* processing done. 
                                                # Prevents re-reading from disk.
        self.image_processed    = np.zeros((1)) # Current image, after all processing is done (bg, scaling, etc)
        self.image_bg_raw       = np.zeros((1)) # The 'raw' background image. No processing done to it. 
                                                # Assume bg is single image; revisit as needed.
        
        self.file_backplane_shortname  = ''     # Shortname for the currently loaded backplane.
        self.file_objectlist_shortname = ''     # Shortname for the currently loaded objectlist.
								
        self.legend             = False         # Store pointer to plot legend here, so it can be deleted
								
        (dim, radii)            = sp.bodvrd('JUPITER', 'RADII', 3)
        self.rj                 = radii[0]      # Jupiter radius, polar, in km. Usually 71492.
        self.is_unwrapped       = False         # Flag: Did we successfully unwrap this image? 
        self.is_loaded          = False         # Flag: Is there a good image loaded?
        
# Now define and startup the widgets
                
        self.path_settings = "/Users/throop/python/salt_interact_settings/"

# Create the sliders, for navigation offset
        
        self.slider_offset_dx  = ttk.Scale(master, from_=-450, to=450, orient=tkinter.HORIZONTAL, 
                                   command=self.set_offset) # Slider dx offset
        self.slider_offset_dy  = ttk.Scale(master, from_=-450, to=450, orient=tkinter.VERTICAL, 
                                   command=self.set_offset) # Slider dy offset

# Define labels

        self.var_label_info = tkinter.StringVar()
        self.var_label_info.set("init")
        self.label_info = ttk.Label(master, textvariable = self.var_label_info)
        
        self.var_label_status_io = tkinter.StringVar()
        self.var_label_status_io.set('IO STATUS')
        self.label_status_io = ttk.Label(master, textvariable = self.var_label_status_io) # For the 'save status' text
                   
# Create the buttons
   
        self.button_prev =     ttk.Button(master, text = 'Prev',     command=self.select_image_prev)
        self.button_next =     ttk.Button(master, text = 'Next',     command=self.select_image_next)
        self.button_plot   =   ttk.Button(master, text = 'Plot',     command=self.handle_plot_image)
        self.button_extract =  ttk.Button(master, text = 'Extract',  command=self.handle_extract_profiles)
        self.button_recenter = ttk.Button(master, text = 'Recenter', command=self.handle_recenter)
        self.button_repoint =  ttk.Button(master, text = 'Repoint',  command=self.repoint)
        self.button_ds9 =      ttk.Button(master, text = 'Open in DS9', command=self.open_external)
        self.button_copy =     ttk.Button(master, text = 'Copy name to clipboard', command=self.copy_name_to_clipboard)
        self.button_quit =     ttk.Button(master, text = 'Quit',     command = self.quit)
        self.button_load =     ttk.Button(master, text = 'Load',     command = self.load)
        self.button_save =     ttk.Button(master, text = 'Save',     command = self.save_file)
        self.button_export =   ttk.Button(master, text = 'Export Analysis', command = self.export_analysis)
        
        self.var_checkbutton_extract = tkinter.IntVar()
        self.var_checkbutton_extract.set(self.do_autoextract)        # Values are in fact 1,0  not  True,False. Ugh. 
        self.checkbutton_extract = ttk.Checkbutton(master, text = 'Auto-Extract + Export', 
                                                   command = self.handle_checkbutton, 
                                                   var = self.var_checkbutton_extract)

        
# Create controls for the background subtraction
# NB: Sometimes the current value for OptionMenu is left as blank. If this happens, need to restart ipython console.
      
        self.var_option_bg = tkinter.StringVar()
        self.var_option_bg.set("Polynomial Order:")
        
        self.label_bg =         ttk.Label(master, text='Background Method')
        self.option_bg =        ttk.OptionMenu(master, self.var_option_bg,\
                                    option_bg_default,  # First argument is the default
                                    "Polynomial", "Previous", "Next", "Median Range", "None", "Grp Num Frac Pow", \
                                     "String", command = self.select_bg)

        self.entry_bg=          ttk.Entry(master, width=16)
        self.entry_bg.insert(0, entry_bg_default) # Set the value (e.g., order = "4")
              
# Create the Entry boxes (ie, single-line text)
# Note that for .insert(num, text), 'num' is in format line.character. Line starts at 1; char starts at 0.
# For both .text and .entry .
        
        self.entry_comment = ttk.Entry(master)
        self.entry_comment.insert(0, 'Comment')

# Create the text box to stuff the header into
# Since we are just displaying text, I'm not sure if this is the right widget to use
# http://www.tkdocs.com/tutorial/text.html
        
        self.text_header = tkinter.Text(master, height=15)

        self.text_header.insert(1.0,'HEADER goes here')

# Add the scrollbars for the header data

        self.scrollbar_text_header_y = ttk.Scrollbar(master, orient='vertical', command=self.text_header.yview)
        self.text_header['yscrollcommand'] = self.scrollbar_text_header_y.set
        
# Create the Listbox for Groups (ie, 'Desc') and load it up
# In theory I should use listvar and StringVar() for this, but they seem to not work for arrays with 
# spaces in them! So I have to load up the elments individually.
# Need to do exportselection=False. Otherwise, only one listbox will be 'active' at a time, and won't be
# able to see value of current setting in both. 
#   http://stackoverflow.com/questions/756662/using-multiple-listboxes-in-python-tkinter

        self.lbox_groups = tkinter.Listbox(master, height=5, selectmode = "browse",
                                           width=25, exportselection=False)
                                           
        for ii,item in enumerate(self.groups.data):    # "groups.data" gives the name of the group
            self.lbox_groups.insert("end", repr(ii) + ". " + item)

        self.lbox_groups.bind('<<ListboxSelect>>', self.select_group) # Define the event handler for 
# Listbox is tied to an event handler, *not* to the normal command= thing.

# Create the listbox for files

        files_short = self.t_group['Shortname']
        num_files = np.size(files_short)

        self.lbox_files = tkinter.Listbox(master, height=20, selectmode = "browse",
                                            width = 30, exportselection=False)

# Populate the listbox for files

        for index_image in range(np.size(files_short)):    
            line_new = self.make_info_line(self.index_group, index_image)
            print("Adding line: {}".format(line_new))
            self.lbox_files.insert("end", line_new)
                  
        self.lbox_files.bind('<<ListboxSelect>>', self.select_image) # Define event handler

# Call select_group, to load the contents of the current group (from index_group) into the file list

#        self.select_group(self)
        
# Add the scrollbars for the file listbox

        self.scrollbar_files_y = ttk.Scrollbar(master, orient='vertical', command=self.lbox_files.yview)
        self.lbox_files['yscrollcommand'] = self.scrollbar_files_y.set
        
# Set up Plot 1 : Main image window
        
# add_subplot: " kwargs are legal Axes kwargs. The Axes instance will be returned."
# http://matplotlib.org/api/figure_api.html
# http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
     
        self.fig1 = Figure(figsize = (7,7))    # <- this is in dx, dy... which is opposite from array order!
        junk = self.fig1.set_facecolor(self.bgcolor)

        self.ax1 = self.fig1.add_subplot(1,1,1, 
                                    xlabel = 'X', ylabel = 'Y', 
                                    label = 'Image',
                                    xlim = (0,1023),
                                    ylim = (1023,0)) # Return the axes
 
        self.canvas1 = FigureCanvasTkAgg(self.fig1,master=master)
        self.canvas1.show()
        self.ax1.imshow(hbt.read_lorri('')) # Put up an empty frame, if file = ''
        
# Set up Plot 2 : Radial / Azimuthal profiles
        
        self.fig2 = Figure(figsize = (7,3))
        junk = self.fig2.set_facecolor(self.bgcolor)

        self.ax2 = self.fig2.add_subplot(1,1,1, 
                                    xlabel = 'Radius or Azimuth', ylabel = 'Intensity', 
                                    label = 'Plot') # Return the axes

        self.canvas2 = FigureCanvasTkAgg(self.fig2,master=master)
        self.canvas2.show()  

# Set up Plot 3 : Unwrapped ring
        
        self.fig3 = Figure(figsize = (7,4))
        junk = self.fig3.set_facecolor(self.bgcolor)

        self.ax3 = self.fig3.add_subplot(1,1,1, 
                                    xlabel = 'Azimuth', ylabel = 'Intensity', 
                                    label = 'Plot') # Return the axes
 
        self.canvas3 = FigureCanvasTkAgg(self.fig3,master=master)
        self.ax3.imshow(hbt.read_lorri(''))
        self.canvas3.show()  

# Set up Plot 4 : Debugging / Misc window
     
        self.fig4 = Figure(figsize = (7,3))
        junk = self.fig4.set_facecolor(self.bgcolor)

        self.ax4 = self.fig4.add_subplot(1,1,1, 
                                    xlabel = 'X', ylabel = 'Y', 
                                    label = 'Image') # Return the axes
 
        self.canvas4 = FigureCanvasTkAgg(self.fig4,master=master)
        self.canvas4.show()
        plot4 = self.ax4.plot([1,2], [3,5]) # Put up an empty frame, if file = ''
        
# Put objects into appropriate grid positions

# Column 1 (LHS) -- numbered as column 1 .. 2

        self.lbox_groups.grid(  row=1, column=1, sticky = 'news', rowspan=1, padx=10)
        self.lbox_files.grid(   row=2, column=1, sticky = 'news',rowspan=2, padx=10)
        self.entry_comment.grid(row=4, column=1, sticky = 'new')
        self.label_info.grid(   row=5, column=1, sticky='new', padx=10, )
        self.text_header.grid(  row=6, column=1, rowspan=6, sticky = 'new', padx=10, pady=10)

        self.scrollbar_files_y.grid(       row=2, column=2, rowspan=2, sticky = 'ns')
        self.scrollbar_text_header_y.grid( row=6, column=2, rowspan=6, sticky = 'ns')
               
# Column 2 (Center) -- numbered as column 10 .. 14

        self.slider_offset_dy.       grid(row=1,  column=10, rowspan=3, sticky='ns') 

        self.canvas1.get_tk_widget().grid(row=1,  column=11, columnspan=4, rowspan=3, sticky = 'nsew')
        self.slider_offset_dx.       grid(row=5,  column=11, columnspan=4, sticky='new')

        self.label_bg.               grid(row=6,  column=11)
        self.option_bg.              grid(row=6,  column=13, sticky='we')
        self.entry_bg.               grid(row=6,  column=14)

        self.label_status_io.        grid(row=7,   column=11) 

        self.button_recenter.        grid(row=8,  column=11)
        self.button_plot.            grid(row=8,  column=13)
        self.button_extract.         grid(row=8, column=14)
    
        self.button_prev.            grid(row=9,  column=11, columnspan=1, sticky='we')
        self.button_next.            grid(row=9,  column=13, columnspan=1, sticky='we')
        self.checkbutton_extract.    grid(row=9,  column=14)

        self.button_load.            grid(row=10,  column=11)
        self.button_save.            grid(row=10,  column=13)
        self.button_export.          grid(row=10,  column=14)
        
        self.button_quit.            grid(row=11,  column=11)
        self.button_ds9.             grid(row=11,  column=13)
        self.button_copy.            grid(row=11,  column=14)
 
# Column 3 (RHS) -- numbered as column 20
        
        self.canvas2.get_tk_widget().grid(row=1, column=20, rowspan = 2)
        self.canvas3.get_tk_widget().grid(row=3, column=20, rowspan = 2)
        self.canvas4.get_tk_widget().grid(row=5, column=20, rowspan = 9)
        
# Define some keyboard shortcuts for the GUI
# These functions must be defined as event handlers, meaning they take two arguments (self and event), not just one.

        master.bind('q', self.quit_e)
        
# Finally, now that GUI is arranged, load the first image, and its backplane.
                
        self.refresh_gui()
        self.load_image()
        self.process_image()
        self.plot_image()

        if (self.do_autoextract == 1):
            self.extract_profiles()

# =============================================================================
# Return a line of info about a file (exptime, mode, filename, time, etc.) 
# This goes into the large 'listbox' table to click on.
# =============================================================================
            
    def make_info_line(self, index_group, index_image):  
        
        import humanize
        import os
        import datetime
        
        # The group is passed, because sometimes it will be current group, and sometimes not
        
        # Get the modtime for the output pickle file
        
        file_analysis = self.get_export_analysis_filename(index_group, index_image)
        
        groupmask = self.t['Desc'] == self.groups[index_group]
        t_group = self.t[groupmask]
        
        print("For index {}, file = {}".format(index_image, file_analysis))
        
        if (os.path.isfile(file_analysis)):               # If the analysis file exists
            time_file = os.path.getmtime(file_analysis)
            time_now  = time.time()
            dt        = time_now - time_file
            timestr   = humanize.naturaltime(dt)          # Create a string like '3 minutes ago'
                                                          # This doesn't update automatically after each edit -- 
                                                          # just after a reload of the group. Still, better than
                                                          # nothing.

        else:                                               # If no file found
            timestr   = '--'
            
        s = '{:3}.   {}   {:6.2f}   {}  {}    {}'.format(   # Create the output string.
             repr(index_image), 
             t_group['Shortname'][index_image].
                 replace('_opnav.fit', '').
                 replace('0x633_sci_', '').
                 replace('0x630_sci_','').
                 replace('lor_',''),       # Shorten filename a bit
             t_group['Exptime'][index_image], 
             t_group['Format'][index_image],
             t_group['UTC'][index_image].split('.')[0],           # Remove fractional seconds                 
             timestr)
        return s

# =============================================================================
# Get the name of the analysis file to export
# =============================================================================

    def get_export_analysis_filename(self, index_group = None, index_image = None):

        if (index_image is None):
            index_image = self.index_image  # Use the current image, unless one is passed

        if (index_group is None):
            index_group = self.index_group

#        else:
                           # Use the passed-in image name
        
        groupmask = self.t['Desc'] == self.groups[index_group]
        t_group = self.t[groupmask]  # 
        
        t = t_group[index_image]  # Grab this, read-only, since we use it a lot.

        dir_export = '/Users/throop/data/NH_Jring/out/'
        file_export = dir_export + t['Shortname'].replace('_opnav', '').replace('.fit', '_analysis.pkl')

        return(file_export)
        
#==============================================================================
# Unwrap ring image in ra + dec to new image in lon + radius.
#==============================================================================

    def unwrap_ring_image(self):

        self.diagnostic('unwrap_ring_image()')
        
        (numrad, rj_array) = sp.bodvrd('JUPITER', 'RADII', 3) # 71492 km
        rj = rj_array[0]
        r_ring_inner = 126000   # Follow same limits as in Throop 2004 J-ring paper fig. 7
        r_ring_outer = 132000

        num_bins_azimuth = 300    # 500 is OK. 1000 is too many -- we get bins ~0 pixels
        num_bins_radius  = 300
             
 
# Check if the backplane is loaded already. Load it iff it is not loaded
    
        if (self.file_backplane_shortname != self.t_group['Shortname'][self.index_image]):
            self.load_backplane()

# Grab the already-loaded object mask for the current image and pointing. 

        mask_objects    = self.objectmask

# Unwrap the ring image. The input to this is the processed ring (ie, bg-subtracted)

        dx = self.offset_dx # Navigational offset. Ideally this would be sub-pixel, although np.roll() only allow ints.
        dy = self.offset_dy
        
        try:
            (im_unwrapped, mask_unwrapped, bins_radius, bins_azimuth) = \
                                      nh_jring_unwrap_ring_image(self.image_processed, 
                                                                 num_bins_radius, (r_ring_inner, r_ring_outer),
                                                                 num_bins_azimuth, 
                                                                 self.planes, 
                                                                 dx=-dx, dy=-dy, 
                                                                 mask=mask_objects)
            self.is_unwrapped = True
        
        except ValueError:
#            self.is_unwrapped = False
            print('Cannot unwrap ring image -- no valid points')
            return
            
# Save the variables, and return
        
        self.image_unwrapped    = im_unwrapped     # NB: This has a lot of NaNs in it.
        self.radius_unwrapped   = bins_radius      # The bins which define the unwrapped coordinates
        self.azimuth_unwrapped  = bins_azimuth     # The bins which define the unwrapped coordinates
        self.mask_unwrapped     = mask_unwrapped   # The object mask, unwrapped. 1 = [object here]

        return

#==============================================================================
# Extract ring profiles from the data image. Also plots them.
#==============================================================================
        
    def extract_profiles(self):
        
        self.diagnostic('extract_profiles()')
        
        # Compute additional quantities we need

        (vec, lt)        = sp.spkezr('New Horizons', 
                            self.t_group['ET'][self.index_image], 'IAU_JUPITER', 'LT', 'Jupiter')
        (junk, lon, lat) = sp.reclat(vec[0:3])
        ang_elev         = np.abs(lat)          # Elevation angle (aka 'B')  
        ang_emis         = math.pi/2 - ang_elev     # Emission angle (ie, angle down from normal) 
        mu               = abs(math.cos(ang_emis))  # mu. See definitions of all these Throop 2004 @ 63 

        width_ring  = 6500 # in km. This is a constant to s by -- not the same as actual boundaries.
        
        ew_edge     = [120000, 130000]  # Integrate over this range. This is wider than the official ring width.

        self.unwrap_ring_image()        # Creates arrays self.{image,radius,azimuth}_unwrapped
                                        # This also sets the is_unwrapped flag
                                        
        if (self.is_unwrapped == False):
#            self.canvas3.delete('all') # Doesn't work
#            self.canvas4.delete('all')
            return
        
        self.ang_elev   = ang_elev
        
        # Load the unwrapped image. ** I should rename the local vars to be same as the self.vars! Only historical.
        
        dn_grid       = self.image_unwrapped # This is only set if this image was navigated, and could be extracted.
        bins_radius   = self.radius_unwrapped
        bins_azimuth  = self.azimuth_unwrapped

#==============================================================================
# Extract radial and azimuthal profiles.
#==============================================================================
        
#        We might use a lot of different cases here. Rather than define each as individual variables, do as dict

        profile_azimuth  = {'inner' : np.array([0]), 'core'   : np.array([0]), 'outer' : np.array([0])}
        profile_radius   = {'full'  : np.array([0]), 'center' : np.array([0]), 'core' : np.array([0])}
        
        range_of_radius  = {'inner' : (126500,127500), 'core' : (127500,129500), 'outer' : (130000,131000)} # for az
        range_of_azimuth = {'full'  : 1, 'center' : 0.25, 'core' : 0.1}                                    # for radial 
        
# Make radial profiles
        
        for key in profile_radius:
            
            profile_radius[key] = nh_jring_extract_profile_from_unwrapped(self.image_unwrapped, 
                                                  self.radius_unwrapped,   # Array defining bins of radius
                                                  self.azimuth_unwrapped,  # Array defining bins of azimuth
                                                  range_of_azimuth[key],        # ie, range of az used for rad profile
                                                  'radius',
                                                  mask_unwrapped = self.mask_unwrapped)   

# Make azimuthal profiles

        for key in profile_azimuth:
            
            profile_azimuth[key] = nh_jring_extract_profile_from_unwrapped(self.image_unwrapped, 
                                                  self.radius_unwrapped,   # Array defining bins of radius
                                                  self.azimuth_unwrapped,  # Array defining bins of azimuth
                                                  range_of_radius[key],        # range of rad used for az profile
                                                  'azimuth',
                                                  mask_unwrapped = self.mask_unwrapped)   

# Create a final 'net' profile, by subtracting the inner and outer background levels from the core profile

        profile_azimuth['net'] = profile_azimuth['core'] - (profile_azimuth['inner'] + profile_azimuth['outer'])/2
        range_of_radius['net'] = range_of_radius['core']
        
        plt.rcParams['figure.figsize'] = 16,10
        plt.rcParams['figure.figsize'] = 10,5

# Save the profiles to self.

        self.profile_radius   = profile_radius
        self.profile_azimuth  = profile_azimuth
        self.range_of_radius  = range_of_radius
        self.range_of_azimuth = range_of_azimuth
        
#==============================================================================
# Convert extracted values from DN, into photometric quantities
#==============================================================================


        rsolar = 266400.  #   RSOLAR  =  266400.000000 / Conv to radiance for solar source. From LORRI FITS file.
                          #   XXX As per HAL, we should not trust the FITS header. Instead, list explicitly.
                          
        pixfov = 0.3 * hbt.d2r * 1e6 / 1024  # This keyword is missing. "Plate scale in microrad/pix "

# Convert into I/F

        profile_radius_iof       = profile_radius.copy()
        profile_radius_iof_norm  = profile_radius.copy()
        
        for key in profile_radius:
            
            profile_radius_iof[key] = self.dn2iof(profile_radius[key], self.t_group['Exptime'][self.index_image],\
                                     pixfov, rsolar)

            profile_radius_iof_norm[key] = profile_radius_iof[key] * mu

        ew_edge_bin = hbt.x2bin(ew_edge,bins_radius)
            
        dr = np.roll(bins_radius,-1) - bins_radius # Bin width, in km
        dr[-1] = 0

# Convert into equivalent width
        
#        ioprofile_radius_iof_norm  # Just a shorter alias
#        ew_norm  = np.sum((profile_radius_full_iof_norm * dr)[ew_edge_bin[0] : ew_edge_bin[1]]) 
                                         # Normalized EW from top
#        ew_mean  = ew_norm / width_ring                                           # Mean normalized EW
        
#        taupp    = profile_radius_full_iof_norm * 4 * mu  # Radially averaged tau_omega0_P
         
#==============================================================================
#  Plot the remapped 2D images
#==============================================================================

        extent = [bins_azimuth[0], bins_azimuth[-1], bins_radius[0], bins_radius[-1]]

        f = (np.max(bins_radius) - np.min(bins_radius)) / (np.max(bins_azimuth) - np.min(bins_azimuth))
        aspect = 0.5 * 1/f  # Set the aspect ratio

        self.ax3.clear()        

        self.ax3.imshow(self.stretch(dn_grid), extent=extent, aspect=aspect, origin='lower')

        self.ax3.set_ylim([bins_radius[0], bins_radius[-1]])
        self.ax3.set_xlim([bins_azimuth[0], bins_azimuth[-1]])
        
        self.ax3.set_xlabel('Azimuth [radians]')
        self.ax3.set_ylabel('Radius [$R_J$]')
        
        self.canvas3.show()
            
#==============================================================================
# Plot the radial and azimuthal profiles
#==============================================================================

        alpha_profiles = (0.3, 1) # Opacity to use for curves, depending on which line it is. 
                             # I match with a boolean test: order = (False, True)
        
# Plot azimuthal profile

        self.ax2.clear()  # Clear lines from the current plot. 
        
        dy = 2            # Vertical offset between curves
                          # Set the y limit to go from minimum, to a bit more than 90th %ile
                          # The idea here is to clip off flux from a moon, if it is there.
        
        frac_edge = 0.1   # How much of the edge of these profiles do we ignore, when setting the range?
        
        vals_central = []
           
        for key in profile_azimuth:      
            self.ax2.plot(bins_azimuth, profile_azimuth[key], 
                          label = key + ' ' + repr(range_of_radius[key]), 
                          alpha = alpha_profiles[key == 'net'])

            bin_left  = int(np.size(bins_azimuth)*frac_edge)
            bin_right = int(np.size(bins_azimuth)*(1-frac_edge))
            
            vals_central_i = profile_azimuth[key][bin_left:bin_right]
            vals_central.append(vals_central_i)
            
        self.ax2.set_title('Azimuthal Profile')
        self.ax2.set_xlabel('Azimuth [radians]')
        self.ax2.set_xlim([bins_azimuth[0], bins_azimuth[-1]])
        self.ax2.legend(loc = 'lower left')
        
        # Set the vertical scale. This takes some subtlety to do correctly, since there could be some invalid etc.
        # Also, the values on the extremes of the range (first/last 10%, eg) are most likely to be wrong, since 
        # they are the ones where the polynomial subtraction could go crazy. So, we want to take the minmax from 
        # the inner (say) 80% of the range, using all three curves.
        
#        profile_azimuth_all
        
        self.ax2.set_ylim(hbt.mm(vals_central))

        self.canvas2.show()

        self.ax4.clear()  # Clear lines from the current plot.
        
        plt.rcParams['figure.figsize'] = 16,10

        DO_PLOT_IOF = False # Do we plot I/F, or DN

# Plot radial profile
        
        alpha_profiles = (0.5, 0.5) # Opacity to use for curves, depending on which line it is. 

        linestyle_ring   = 'dotted'             # Ring linestyle
        alpha_ring       = 0.35
        linewidth_ring   = 0.4
        color_ring       = 'blue'
        
        
        if (DO_PLOT_IOF == False):

            for key in profile_radius: 
                self.ax4.plot(bins_radius/1000, profile_radius[key], label = key + ' ' + repr(range_of_azimuth[key]), 
                              alpha = alpha_profiles[key == 'core'])

            self.ax4.set_xlabel(r'Radial Profile      Radius [1000 km]    phase = {:.1f} deg'
                                .format(hbt.r2d * np.mean(self.planes['Phase'])))

            self.ax4.set_xlim([126,130])
            self.ax4.legend(loc='upper left')
            
            # Plot a second axis, in RJ    
            
            ax41 = self.ax4.twiny()
            ax41.set_xlim(list(hbt.mm(bins_radius/1000/71.4)))
            
            # Draw vertical lines for the ring locations
            
            self.ax4.axvline(x=self.a_ring_inner_km/1000, linestyle=linestyle_ring, color=color_ring,# alpha=alpha_ring, 
                             linewidth=linewidth_ring)
            
            self.ax4.axvline(x=self.a_ring_outer_km/1000, linestyle=linestyle_ring, color=color_ring,# alpha=alpha_ring, 
                             linewidth=linewidth_ring)
            

        if (DO_PLOT_IOF == True):     
            
            self.ax4.plot(bins_radius / 1000, taupp, label = 'Method 1')
            self.ax4.set_xlabel('Radius [1000 km]')
            self.ax4.set_ylabel(r'$\tau \varpi_0 P(\alpha)$')
             
            self.ax4.set_xlabel('Radial Profile      Radius [1000 km]')

            self.ax4.set_xlim(list(hbt.mm(bins_radius/1000)))
            self.ax4.legend(loc='upper left')
    
            self.ax4.plot(self.bins_radius, self.profile_radius)
            ax41 = self.ax4.twiny()
            ax41.set_xlim(list(hbt.mm(bins_radius/1000/71.4)))
            
        self.canvas4.show()
        
##########
# Recenter the image. That is, this resets the slider positions to (0,0). 
# Ideally we would have a display and read these out dynamically, but this will do the job.
##########

    """
    Handle the 'Recenter image' button.
    Plots image when done.
    """
				
    def handle_recenter(self):
        
        self.diagnostic("handle_recenter()")
        
        dx = 0
        dy = 0

        # Set the slider position
        
        self.slider_offset_dx.set(dx)
        self.slider_offset_dy.set(dy)
        
        # Set the internal position
        
        self.t_group['dx_offset'][self.index_image] = dx
        self.t_group['dy_offset'][self.index_image] = dy
        
        # Redraw objects
        
        self.clear_objects()
        self.plot_objects()
        
        return True								
							
#########
# Navigate the image
#########

    def navigate(self):        
                    
        return True

##########
# Save GUI settings. 
##########
    """
    Take all settings from the GUI, and put them into the proper table variables
    This does *not* save to disk. It saves the Tk widget settings into variables, 
    which can then be saved later.
    """
 
    def save_gui(self):

#        print("save_gui()")
        
# Save the current values
         
        comment = self.entry_comment.get()
        self.t_group['Comment'][self.index_image] = comment

        # Get the slider positions
        
        self.t_group['dx_offset'][self.index_image] = self.slider_offset_dx.get()
        self.t_group['dy_offset'][self.index_image] = self.slider_offset_dy.get()

        # Get the background subtraction settings

        self.t_group['bg_method'][self.index_image] = self.var_option_bg.get()
        self.t_group['bg_argument'][self.index_image] = self.entry_bg.get()
        
# Copy the entire current 'group' back to the main table

        self.t[self.groupmask] = self.t_group
        
##########
# Change Image, and go to a new one
##########
#
# index_image     : old image has the current image (within the current group)
# index_image_new : new image
                                # has the new image we want to show (within the current group)
# This does not change the group -- just the image      

    def change_image(self, refresh=True):
   
#        print("change_image")

        # Save current GUI settings into variables
        
        self.save_gui()
        
# Update the index, to show that new image is active
        
        self.index_image = self.index_image_new
        self.index_group = self.index_group_new

#  Update the frequently used 't_group' list to reflect new group

        name_group = self.groups[self.index_group]
        self.groupmask = self.t['Desc'] == name_group        
        self.t_group = self.t[self.groupmask]
        self.num_images_group = np.size(self.t_group)							
								
# Redraw the GUI for the new settings

        self.refresh_gui()

# Load, process, and display the image
								
        self.load_image()
					
        self.process_image()
								
        self.plot_image()
                
# Extract radial / azimuthal profiles, and export results, if desired

        if (self.do_autoextract == 1):
            self.extract_profiles() 
            self.export_analysis()       # Write results of new image to disk. Otherwise, it is easy to forget.

#==============================================================================
# Load new image from disk
#==============================================================================

    """
    Load current image from disk
    """
    
    def load_image(self):

        self.diagnostic("load_image()")

# autozoom: if set, and we are loading a 4x4 image, then scale it up to a full 1x1.

        file = self.t_group[self.index_image]['Filename']  

        self.image_raw = hbt.read_lorri(file, frac_clip = 1., bg_method = 'None', autozoom=False)
        self.diagnostic("Loaded image: " + file)

# Load info from header

        header = hbt.get_image_header(file)
        self.et = header['SPCSCET']
        self.utc = sp.et2utc(self.et, 'C', 0)
        
# Load the backplane as well. We need it to flag satellite locations during the extraction.
        
        DO_LOAD_BACKPLANE = True
        if DO_LOAD_BACKPLANE:
            self.load_backplane()
            
# Load the 'objectlist' in this frame as well. This is the precomputed list of stars and satellites.
# It is used as a mask.

        DO_LOAD_OBJECTLIST = True
        if DO_LOAD_OBJECTLIST:
            self.load_objectlist()
        
# Set a flag to indicate that a valid image has been loaded
            
        self.is_loaded = True
        
#==============================================================================
# Change Group
#==============================================================================

    def select_group(self, event):

        self.diagnostic("select_group")

        index_group = (self.lbox_groups.curselection())[0]
        name_group = self.groups[index_group]

# Create the listbox for files, and load it up. This is the graphic 'data table' in the GUI.
# NB: We are using a local copy of t_group here -- not self.t_group, 
# which still reflects the old group, not new one.

        groupmask = self.t['Desc'] == name_group        
        t_group = self.t[groupmask]
        num_images_group = np.size(t_group)
        print("Group #{} = {} of size {} hit!".format(index_group, name_group, self.num_images_group))

# Now look up the new group name. Extract those elements.        

        files_short = t_group['Shortname']

        self.lbox_files.delete(0, "end")
 
# For each line in the GUI table, create it, and then load it.

        for index_image in range(np.size(files_short)):    
            line_new = self.make_info_line(index_group, index_image)
#            print("Adding line: {}".format(line_new))
            self.lbox_files.insert("end", line_new)
                 
# Then set the group number, and refresh screen.

        self.index_image_new = 0      # Set image number to zero
        self.index_group_new = index_group  # Set the new group number
#        print(" ** Calling change")
        self.change_image()

##########
# Select new image, based on user click
##########
# This is called when user clicks on a new image number in list.

    def select_image(self, event):

        self.diagnostic("select_image", False)
        
        index = (self.lbox_files.curselection())[0]
        name = self.t_group['Shortname'][index]
        print()
        self.diagnostic("selected = " + repr(index) + ' ' + name)
        self.index_image_new = index
        self.index_group_new = self.index_group  # Keep the same group
        self.change_image()

##########
# Refresh the status bar
##########

    def refresh_statusbar(self):

# Set the Label widget value with a useful string
# Note that this is sometimes just not printed at all.  This is usually caused by a non-closed Tk window.
# To fix, quit IPython console from within spyder. Then restart.
								
        num_in_group = np.size(self.t_group)
        num_groups = np.size(self.groups)
                
        s = 'Group ' + repr(self.index_group) + '/' + repr(num_groups) + \
            ', Image ' + repr(self.index_image) + '/' + repr(num_in_group) + ': ' + \
          self.t_group['Shortname'][self.index_image] + '\n'
        
        if (self.t_group['is_navigated'][self.index_image]):
            s += ('Nav, ' + \
            'd{xy}_opnav = (' + repr(self.t_group['dx_opnav'][self.index_image]) + ', ' + \
                                repr(self.t_group['dy_opnav'][self.index_image]) + '), ' + \
            'd{xy}_user =  (' + repr(int(self.slider_offset_dx.get())) + ', ' \
                              + repr(int(self.slider_offset_dy.get())) + ')')
            
        self.var_label_info.set(s)

        return 0
        
##########
# Update UI for new settings.
##########
        """
A general-purpose routine that is sometimes called after (e.g.) switching images.
the image number. This updates all of the GUI elements to reflect the new image, based on 
the internal state which is already correct. This does *not* refresh the image itself.
        """

    def refresh_gui(self):
         
# Load and display the proper header

        filename = self.t_group['Filename'][self.index_image]
        
        header = hbt.get_image_header(filename, single=True)
        self.text_header.delete(1.0, "end")
        h2 = re.sub('\s+\n', '\n', header)  # Remove extra whitespace -- wraps around and looks bad
        self.text_header.insert(1.0, h2)

# Make sure the proper comment is displayed. (It won't be if we have just changed images.)

        self.entry_comment.delete(0, 'end') # Clear existing comment
        self.entry_comment.insert(0, self.t_group['Comment'][self.index_image])
        
# Set the proper Listbox 'Group' settings. Make sure the current item is visible, and selected.
# It's pretty easy for the Listbox to get confused, and highlight two-at-once, or underline one and highlight
# another, etc. This just cleans it up. Maybe there are settings to change here too.
    
        self.lbox_groups.see(self.index_group)
        self.lbox_groups.selection_clear(0, 'end')
        self.lbox_groups.selection_set(self.index_group) # This seems to not work always
        self.lbox_groups.activate(self.index_group)

# Set the proper Listbox 'image' settings.

        self.lbox_files.see(self.index_image)
        self.lbox_files.selection_clear(0, 'end')
        self.lbox_files.selection_set(self.index_image)
        self.lbox_files.activate(self.index_image)

# Set the slider positions

        self.slider_offset_dx.set(self.t_group['dx_offset'][self.index_image])
        self.slider_offset_dy.set(self.t_group['dy_offset'][self.index_image])

#        self.slider_offset_dy.set(self.t_group['dy_offset'][self.index_image_new]) # Set dy

# Set the background subtraction settings

        self.var_option_bg.set(self.t_group['bg_method'][self.index_image])
        self.entry_bg.delete(0, 'end') # Clear existing value
        self.entry_bg.insert(0, self.t_group['bg_argument'][self.index_image]) # Set the value (e.g., order = "4")
        
# Set the statusbar

        self.refresh_statusbar()
        
        return 0

#==============================================================================
# Extract profiles - button handler
#==============================================================================

    def handle_extract_profiles(self):
        """
        Extract and plot radial / azimuthal profiles. Because plotting params
        have likely changed, we replot as well.
        Identical to handle_plot_image, but forces the extraction.								
        """
        
        self.save_gui()
        self.process_image()
        self.plot_image()

        self.extract_profiles()

#==============================================================================
# Apply image processing to image. Stray light, polynomial subtraction, etc.
#==============================================================================

    def process_image(self):
        
        self.diagnostic("process_image()")
        
        method = self.var_option_bg.get()
        argument = self.entry_bg.get()
        
        self.image_processed = hbt.nh_jring_process_image(self.image_raw, \
                                                          method, argument, 
                                                          index_group = self.index_group, 
                                                          index_image = self.index_image)

        # Remove cosmic rays
        
        self.image_processed = hbt.decosmic(self.image_processed)
        
#==============================================================================
# Plot image - button handler
#==============================================================================

    def handle_plot_image(self):
        """
        Plot the image. Reprocess to reflect any changes from GUI.
        Plots objects and extracts profiles, as available.								
        """
        
        self.save_gui()
        self.process_image()
        self.plot_image()

        if (self.do_autoextract):
            self.extract_profiles()
            self.export_analysis()
							

#==============================================================================
# Plot image
#==============================================================================

    def plot_image(self):
        """
        Plot the image. It has already been loaded and processed.
        This is an internal routine, which does not process.	
        This plots the image itself, *and* the objects, if they exist.							
        """
        
        self.diagnostic("plot_image()")
        
        # Clear the image

        self.ax1.clear()
        
        # Set the color map
        
        plt.rc('image', cmap='Greys_r')
        
# Plot the image itself.
        
        # Render the main LORRI frame
        # *** In order to get things to plot in main plot window, 
        # use self.ax1.<command>, not plt.<command>
        # ax1 is an instance of Axes, and I can call it with most other methods (legend(), imshow(), plot(), etc)
                       
        self.ax1.imshow(self.stretch(self.image_processed))
        
        # Disable the tickmarks from plotting

        self.ax1.get_xaxis().set_visible(False)
        self.ax1.get_yaxis().set_visible(False)

        # Set image size (so off-edge stars are clipped, rather than plot resizing)

        # Display it and scale spatially, whether it is 1x1 or 4x4
        
        self.ax1.set_xlim([0,hbt.sizex(self.image_processed)-1])  # This is an array and not a tuple. 
                                                                  # Beats me, like so many things with mpl.
        self.ax1.set_ylim([hbt.sizex(self.image_processed)-1,0])
          
        # Draw the figure on the Tk canvas
        
        self.fig1.tight_layout() # Remove all the extra whitespace -- nice!
        self.canvas1.draw()
                        
        self.plot_objects()
        
        return 0

##########
# Clear all objects from plot
##########
        # Remove all of the current 'lines' (aka points) from the plot. This leaves the axis and the image
        # preserved, but just removes the lines. Awesome. Using ax1.cla() will clear entire 'axis', including image.
        # Q: Does this remove legend? Not sure.

    def clear_objects(self):
    		
            lines = self.ax1.get_lines()
            for line in lines:
                line.remove()		# Not vectorized -- have to do it one-by-one


##########
# Plot Objects
##########

    def plot_objects(self):
        """
           Plot rings and stars, etc, if we have them.
      """

        self.diagnostic("plot_objects()")
        DO_PLOT_STARS      = True
        DO_PLOT_SATS       = True
        DO_PLOT_RING_INNER = True
        DO_PLOT_RING_OUTER = True

        color_stars_phot = 'red'            # Color for stars found photometrically
        color_stars_cat  = 'lightgreen'     # Color for stars in catalog  
        alpha_stars_cat  = 0.5
        color_sats       = 'pink'
        alpha_sats       = 0.5
        
        marker_ring      = 'None'           # Must be capitalized
        linestyle_ring   = 'dotted'             # Ring linestyle
        alpha_ring       = 0.35
        linewidth_ring   = 0.4
        color_ring       = 'blue'
        
        abcorr           = 'LT+S'           # I have been using LT here. But no reason not to use LT+S I think?
        
        num_pts_ring     = 300 
            
        # If we have them, draw the stars on top

        t = self.t_group[self.index_image]  # Grab this, read-only, since we use it a lot.
                                            # We can reference table['col'][n] or table[n]['col'] - either OK

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = WCS(t['Filename'])                  # Look up the WCS coordinates for this frame
         
        filename = self.t_group['Filename'][self.index_image]
#        filename = str(filename)  # What was this supposed to do?
        							          
# Get the user offset position

        dx_pix = 0
        dy_pix = 0
            
        dx = t['dx_opnav'] + self.slider_offset_dx.get()
        dy = t['dy_opnav'] + self.slider_offset_dy.get()

# Plot the stars -- catalog 

        stars = self.objectlist[(self.objectlist['name'] == 'star')]
            
        self.ax1.plot(stars['x_pix'] + dx, stars['y_pix'] + dy, 
                     marker='o', ls='None', 
                     color=color_stars_cat, alpha = alpha_stars_cat, ms=12, mew=1, label = 'Stars, Catalog')
                     
# Get position of satellites and plot them

        sats = self.objectlist[ self.objectlist['name'] != 'star' ]
            
        for i in range(np.size(sats)):
                
# Plot the sats                
            self.ax1.plot(sats['x_pix'][i] + dx, sats['y_pix'][i] + dy, 
                     marker='o', ls='None', 
                     color=color_sats, alpha = alpha_sats, ms=12, mew=1)
# Label each one
                
            self.ax1.text(sats['x_pix'][i] + dx, sats['y_pix'][i] + dy,
                              sats['name'][i], clip_on = True)
            
# Plot the ring

# Get the ring pixel locations
            
            if (DO_PLOT_RING_OUTER):
                
                x_ring2_pix, y_ring2_pix = hbt.get_pos_ring(self.et, num_pts = num_pts_ring,
                                                name_body='Jupiter', 
                                                radius=self.a_ring_outer_km, units='pixels', abcorr=abcorr, wcs=w)

                self.ax1.plot(x_ring2_pix+dx, y_ring2_pix+dy, marker=marker_ring, color = color_ring, 
                                                ls = linestyle_ring,
                                                linewidth=linewidth_ring, alpha = alpha_ring)

#                self.ax1.plot(x_pos_ring2 + dx, y_pos_ring2 + dy, marker='o', color = 'lightblue', ls = '--',
#                              label = 'Ring, OpNav+User')
                    
            if (DO_PLOT_RING_INNER):

                x_ring1_pix, y_ring1_pix = hbt.get_pos_ring(self.et, num_pts = num_pts_ring, 
                                                name_body='Jupiter', 
                                                radius=self.a_ring_inner_km, units='pixels', abcorr=abcorr, wcs=w)

                self.ax1.plot(x_ring1_pix+dx, y_ring1_pix+dy, marker=marker_ring, color=color_ring, 
                                                ls = linestyle_ring, 
                                                linewidth=linewidth_ring, alpha = alpha_ring, ms=8)

            self.legend = self.ax1.legend()  # Draw legend. Might be irrel since remove() might keep it; not sure.

            self.canvas1.draw()

#==============================================================================
# Generate a mask of all the stars in field
#==============================================================================
            
    def get_mask_objects(self):
        """
        Returns a boolean image, set to True for pixels close to a satellite or star.
        """
        
        # Not sure if this works yet
        
        mask = nh_jring_mask_from_objectlist(self.file_objectlist)
        
        return mask
        
##########
# Load backplane
##########
# Returns True if loaded successfully; False if not
        
    def load_backplane(self):

        self.diagnostic("load_backplane()")
        
        dir_backplanes = '/Users/throop/data/NH_Jring/out/'
        file_backplane = dir_backplanes + self.t_group['Shortname'][self.index_image].replace('.fit', '_planes.pkl')

        # Save the shortname associated with the current backplane. 
        # That lets us verify if the backplane for current image is indeed loaded.

        self.file_backplane_shortname = self.t_group['Shortname'][self.index_image]
				
        if (os.path.isfile(file_backplane)):
            lun = open(file_backplane, 'rb')
            planes = pickle.load(lun)
            self.planes = planes # This will load self.t
            lun.close()
            return True

        else:
            print("Can't load backplane " + file_backplane)
            return False

##########
# Load objectlist
##########
# Return True if loaded successfully; False if not

    def load_objectlist(self):
        
        dir_objectlists = '/Users/throop/data/NH_Jring/out/'
        
        file_objectlist = dir_objectlists + self.t_group['Shortname'][self.index_image].\
            replace('.fit', '_objects.txt')

        # Save the shortname associated with the current objectlist.
        # That lets us verify if the objectlist for current image is indeed loaded.

        self.file_objectlist = file_objectlist
        
        self.file_objectlist_shortname = self.t_group['Shortname'][self.index_image]
				
        if (os.path.isfile(file_objectlist)):
            t = Table.read(file_objectlist, format='csv')
            self.objectlist = t  # This will load self.t

            # While we are at it, also create the object mask, which is a boolean image
        
            self.objectmask = nh_jring_mask_from_objectlist(file_objectlist)
            
            return True

        else:
            print("Can't load objectlist " + file_objectlist)
            return False    
        
##########
# Repoint an image
##########
# Placeholder -- not sure what it will do

    def repoint(self):

        return 0

##########
# Open a file in external viewer such as DS9
##########

    def open_external(self):

        app = '/Applications/SAOImage DS9.app/Contents/MacOS/ds9'
        file = self.t_group['Filename'][self.index_image]

        subprocess.call(['/usr/bin/open', app, file]) # Non-blocking (but puts up a terminal window)
#        subprocess.call([app, file]) # Blocking
        
        return 0

##########
# Choose the background model
##########

    def select_bg(self, event):
        
#        print("bg method = " + self.var_option_bg.get())

        # Grab the GUI settings for the background, and save them into proper variables.

        self.t_group['bg_method'][self.index_group] = self.var_option_bg.get()
        self.t_group['bg_argument'][self.index_group] = self.entry_bg.get()

##########
# Copy current filename (short) to clipboard
##########

    def copy_name_to_clipboard(self):

        str = self.t_group['Shortname'][self.index_image]
        hbt.write_to_clipboard(str)
        print("{} copied to clipboard".format(str))
        
        return 0
        
##########
# Close app and quit
##########
# Define this with one argument if we press a button, or two if we trigger an event like a keypress
# Get error if called with more than it can handle. So we define it here with two.

    def quit_e(self, event):
        self.quit()        
    
    def quit(self):
#        answer = ttk.messagebox.askquestion
        root.destroy() # Q: How does root. get automatically imported here? Not sure.
        root.quit()

##########
# Handle change to a checkbox
##########

    def handle_checkbutton(self):
        
        print("Button = " + repr(self.var_checkbutton_extract.get()))
        self.do_autoextract = self.var_checkbutton_extract.get() # Returns 1 or 0
        
##########
# Load settings
##########
# Load them from the pre-set filename

    def load(self, verbose=True):
        
        lun = open(self.dir_out + self.filename_save, 'rb')
        t = pickle.load(lun)
        self.t = t # All of the variables we need are in the 't' table.
        lun.close()

        if (verbose):
            tkMessageBox.showinfo("Loaded", "File " + self.filename_save + " loaded.")

##########
# Save everything to disk
##########

    def save_file(self, verbose=True):

#    """ 
#    Save all settings to the pre-set filename 
#    """
        # First save the current GUI settings to local variables in case they've not beeen saved yet.
  
        self.save_gui()

        # Put up a message for the user. A dialog is too intrusive.

        self.var_label_status_io.set('SAVING...') # This one doens't work for some reason...
    
        # Write one variable to a file    
    
        print("Saving to {}".format(self.dir_out + self.filename_save))
        
        lun = open(self.dir_out + self.filename_save, 'wb')
        t = self.t
        pickle.dump(t, lun)
        lun.close()
        
        self.var_label_status_io.set('')  # This one works

##########
# Export results for this image to a file
##########

    def export_analysis(self, verbose=True):

        t = self.t_group[self.index_image]  # Grab this, read-only, since we use it a lot.
        dir_export = '/Users/throop/data/NH_Jring/out/'
        
        file_export = self.get_export_analysis_filename() # With no index passed, take current
        
# Prepare the variable to stuff into the single tuple we pass to write into pickle file.
# This is simple, but it means we much be careful to unpack in same order as they were packed.
# Perhaps it's better to write as a dictionary, where we can then extract by name, rather than by order?
        
        image_unwrapped    = self.image_unwrapped
        mask_unwrapped     = self.mask_unwrapped
        radius             = self.radius_unwrapped
        azimuth            = self.azimuth_unwrapped
        profile_radius_dn  = self.profile_radius
        profile_azimuth_dn = self.profile_azimuth
        ang_elev           = self.ang_elev
        exptime            = t['Exptime']
        ang_phase          = np.mean(self.planes['Phase'])
        index_image        = self.index_image
        index_group        = self.index_group
        bg_method          = t['bg_method']
        bg_argument        = t['bg_argument'] 
        et                 = t['ET']
        
        vals = (image_unwrapped,     # Unwrapped image itself
                mask_unwrapped,      # Boolean mask
                radius,              # Axis values for the unwrapped image
                azimuth,             # Axis values for the unwrapped image
                profile_radius_dn,   # Radial profile (several, in a dictionary)
                profile_azimuth_dn,  # Az profile (several, in a dictionary)
                self.range_of_azimuth,
                self.range_of_radius,
                exptime,             # Exposure time
                et,                  # ET
                ang_elev,            # Elevation angle above ring
                ang_phase,           # Phase angle (mean to rings -- not to planet center)
                bg_method,           # Method of background removal
                bg_argument,         # Details of background removal
                index_image,         # Index of image
                index_group)         # Index of image group
        
        lun = open(file_export, 'wb')
        pickle.dump(vals, lun)
        lun.close()
        
        print("Wrote results: {}".format(file_export))
        
##########
# (p)revious image
##########
      
    def select_image_prev(self):
        
        self.index_image_new = (self.index_image - 1) % self.num_images_group
        self.index_group_new = self.index_group
								
        self.change_image()

##########
# (n)ext image
##########
    
    def select_image_next(self):
          
        self.index_image_new = (self.index_image + 1) % self.num_images_group
        self.index_group_new = self.index_group
        self.change_image()

##########
# Process slider positions
##########

    def set_offset(self, event): # Because this is an event, it it passed two arguments
        
# Get the slider positions, for the dx and dy nav offset positions, and put them into a variable we can use

        self.offset_dx = self.slider_offset_dx.get() #
        self.offset_dy = self.slider_offset_dy.get() #
        
#        print("Offset = {}, {}".format(self.offset_dx, self.offset_dy))
        
        # Now undraw and redraw the ring, if it is loaded
        
        if self.is_loaded:

            self.clear_objects()
    
            self.plot_objects()
        
        else:
            print("Not loaded -- no plotting")
            
        return 0

#==============================================================================
# Convert from DN to I/F
#==============================================================================
        
    def dn2iof(self, dn, exptime, pixfov, rsolar):
        
    # Convert DN to I/F. 
    # This is not a general routine -- it's specific to LORRI 1x1 @ Jupiter.
    
        # Calculate LORRI pixel size, in sr

        pixfov = 0.3 * hbt.d2r * 1e6 / 1024  # This keyword is missing. "Plate scale in microrad/pix "

        sr_pix       = (pixfov*(1e-6))**2  # Angular size of each MVIC pixel, in sr
    
        # Look up the solar flux at 1 AU, using the solar spectral irradiance at 
        #    http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
        # or http://www.pas.rochester.edu/~emamajek/AST453/AST453_stellarprops.pdf
    
        f_solar_1au_si     = 1.76                 # W/m2/nm. At 600 nm. 176 is Hal's recommended number. 
        f_solar_1au_cgs    = f_solar_1au_si * 100 # Convert from W/m2/nm to erg/sec/cm2/Angstrom
    
        f_solar_1au        = f_solar_1au_cgs
    
        # Calculate the solar flux at Pluto's distance. [Actual distance is 33.8 AU]
    
        f_solar_jup       = f_solar_1au / (5**2)  # Use 1/r^2, not 1/(4pi r^2)
    
        # Use constants (from Level-2 files) to convert from DN, to intensity at detector.
    
        # This is the equation in ICD @ 62.
    
        i_per_sr    = dn / exptime / rsolar # Convert into erg/cm2/s/sr/Angstrom.
    
        # Because the ring is spatially extended, it fills the pixel, so we mult by pixel size
        # to get the full irradiance on the pixel.
        # *** But, somehow, this is not working. I get the right answer only if I ignore this 'per sr' factor.
    
        DO_OVERRIDE = True  # If true, ignore the missing factor of 'per sr' that seems to be in the conversion
    
        i           = i_per_sr * sr_pix     # Convert into erg/cm2/s/Angstrom
    
        if (DO_OVERRIDE):
            i = i_per_sr
    
        # Calculate "I/F". This is not simply the ratio of I over F, because F in the eq is not actually Flux.
        # "pi F is the incident solar flux density" -- ie, "F = solar flux density / pi"
        # SC93 @ 125
    
        iof = i / (f_solar_jup / math.pi)
        
        return iof

##########
# Diagnostic message
##########

    def diagnostic(self, message, flag=True):
        '''
        Print a diagnostic message. Can turn this off easily, both globally and for individual messages.
        '''
        
        DO_PRINT = False
   
        if (DO_PRINT):
            print(message)
        
        return
    
###########
# Now start the main app
###########

# Start up the widget

root = tkinter.Tk()
app  = App(root)

# set the dimensions and position of the window

root.geometry('%dx%d+%d+%d' % (1610, 850, 2, 2))
root.configure(background='#ECECEC')                # ECECEC is the 'default' background for a lot of ttk widgets, which
                                                    # I figured out using photoshop. Not sure how to change that.

#  window.attributes('-topmost', 1)
#  window.attributes('-topmost', 0)

os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

# Start a profiler, if requested

DO_PROFILE = False

if DO_PROFILE:
    
    file_out_profile = 'profile.bin'
    
    cProfile.run('root.mainloop()', filename=file_out_profile)  # This will call the __init__ function
    
    print("Profiling with cProfile, output to " + file_out_profile)
    
else:    

    while True:
        try:
            root.mainloop()
            break
        except UnicodeDecodeError:  # Fix a TK + py3 bug, which causes Unicode errors on events. 
                    # http://stackoverflow.com/questions/16995969/inertial-scrolling-in-mac-os-x-with-tkinter-and-python
            pass
