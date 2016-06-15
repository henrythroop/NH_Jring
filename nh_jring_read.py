# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

@author: throop

"""

# -*- coding: utf-8 -*-
"""
# Program reads in NH J-ring data and navigates and extracts it.

# Possible features to be added:
#  o Output images as PNG
#  o Select multiple images and display them as an arrayf
#  o Slider to scale the image brightness / contrast
#  o Subtract median? Subract neighbor image?

# Widget program created 18-May-2016 based on previous text version (same name), and salt_interact.py .

@author: throop
"""

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
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import wcsaxes
import hbt
import time

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

#pdb.set_trace()

# First we define any general-purpose functions, which are not part of the class/module.
# We can move these to a different file at some point.


##########
# Find stars in an image
##########
        
def find_stars(im):
    "Locate stars in an array, using DAOphot. Returns N x 2 array with xy positions. No magnitudes."
         
    mean, median, std = sigma_clipped_stats(im, sigma=3.0, iters=5)
    sources = daofind(im - median, fwhm=3.0, threshold=5.*std)
    x_phot = sources['xcentroid']
    y_phot = sources['ycentroid']
        
    points_phot = np.transpose((x_phot, y_phot)) # Create an array N x 2

    return points_phot

##########
# Calc offset between two sets of points
##########

def calc_offset_points(points_1, points_2, shape, plot=False):
    "Calculate the offset between a pair of ordered points -- e.g., an xy list"
    "of star positions, and and xy list of model postns."
    "Returned offset is integer pixels as tuple (dy, dx)."
    
    diam_kernel = 5 # If this is 11, that is too big, and we gt the wrong answer. Very sensitive.

    image_1 = hbt.image_from_list_points(points_1, shape, diam_kernel)
    image_2 = hbt.image_from_list_points(points_2, shape, diam_kernel)
 
    t0,t1 = ird.translation(image_1, image_2) # Return shift, with t0 = (dy, dx). t1 is a flag or quality or something.
    (dy,dx) = t0
    
    if (plot):

        xrange = (0, shape[0]) # Set xlim (aka xrange) s.t. 
        yrange = (shape[1], 0)

        figs = plt.figure()
        ax1 = figs.add_subplot(1,2,1) # nrows, ncols, plotnum. Returns an 'axis'
        ax1.set_aspect('equal') # Need to explicitly set aspect ratio here, otherwise in a multi-plot, it will be rectangular

#        fig1 = plt.imshow(np.log(image_1))
        plt.plot(points_1[:,0], points_1[:,1], marker='o', color='lightgreen', markersize=4, ls='None', label = 'Photometric')
        plt.plot(points_2[:,0], points_2[:,1], marker='o', color='red', markersize=4, ls='None', label = 'Cat')
        plt.legend()
       
        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        plt.title('Raw')
        
        ax2 = figs.add_subplot(1,2,2) # nrows, ncols, plotnum. Returns an 'axis'
        plt.plot(points_1[:,0], points_1[:,1], marker='o', color='lightgreen', markersize=9, ls='None')
        plt.plot(points_2[:,0] + t0[1], points_2[:,1] + t0[0], marker='o', color='red', markersize=4, ls='None')
        ax2.set_aspect('equal')

        plt.xlim(xrange)    # Need to set this explicitly so that points out of image range are clipped
        plt.ylim(yrange)
        plt.title('Shifted, dx=' + repr(dx) + ', dy = ' + repr(dy))
        
        plt.show()
        
    return t0
        
class App:

# Now define the functions that are part of this module (or class, not sure). They are accessed with 
# self.function(), and they have access to all of the variables within self. .

#########
# (N)avigate the image -- plot rings and everything on it
#########
# Things I want to return from here, or at least get set into the global variable table:
#
#  Star positions for the catalog, abcorr. Pixels.
#  Star positions in this frame, abcorr. Pixels.
#  dx and dy offset between these two sets

    def navigate(self):        

        t = self.t_group # We use this a lot, so make it shorter
        d2r = hbt.d2r
        r2d = hbt.r2d
        
        index_image = self.index_image
        
        xlim = (0,1023)
        ylim = (1023,0)
        
# Now look up positions of stars in this field, from a star catalog

        w = WCS(t['Filename'][index_image])                  # Look up the WCS coordinates for this frame
        et = t['ET'][index_image]

        print 'ET[i] =  ' + repr(et)
        print 'UTC[i] = ' + repr(t['UTC'][index_image])
        print 'crval[i] = ' + repr(w.wcs.crval)
#        print 'RA[i] =  ' + repr(t['RAJ2000'][index_image] * d2r) + ' rad'
#        print 'Dec[i] = ' + repr(t['DEJ2000'][index_image] * d2r) + ' rad'
        
        center  = w.wcs.crval  # degrees
        DO_GSC      = True
        DO_USNOA2   = False
        
        if (DO_GSC):
            name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
            ra_stars  = np.array(stars.array['RAJ2000'])*d2r # Convert to radians
            dec_stars = np.array(stars.array['DEJ2000'])*d2r # Convert to radians
            table_stars = Table(stars.array.data)
            
        if (DO_USNOA2):        
            name_cat = 'The USNO-A2.0 Catalogue (Monet+ 1998) 1' # Works but gives stars down to v=17; I want to v=13 
            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
            table_stars = Table(stars.array.data)
            mask = table_stars['Bmag'] < 13
            table_stars_m = table_stars[mask]            

            ra_stars  = table_stars_m['RAJ2000']*d2r # Convert to radians
            dec_stars = table_stars_m['DEJ2000']*d2r # Convert to radians

# Get position of satellites

        name_bodies = np.array(['Metis', 'Adrastea', 'Io'])        
        x_bodies,  y_bodies   = hbt.get_pos_bodies(et, name_bodies, units='pixels', wcs=w)
        ra_bodies, dec_bodies = hbt.get_pos_bodies(et, name_bodies, units='radec', wcs=w)
        
# Get an array of points along the ring

        ra_ring1, dec_ring1 = hbt.get_pos_ring(et, name_body='Jupiter', radius=122000, units='radec', wcs=w)
        ra_ring2, dec_ring2 = hbt.get_pos_ring(et, name_body='Jupiter', radius=129000, units='radec', wcs=w)

                    # Return as radians              
        x_ring1, y_ring1    = w.wcs_world2pix(ra_ring1*r2d, dec_ring1*r2d, 0) # Convert to pixels
        x_ring2, y_ring2    = w.wcs_world2pix(ra_ring2*r2d, dec_ring2*r2d, 0) # Convert to pixels

# Get position of Jupiter, in pixels

#        ra_io, dec_io = get_pos_body(et, name_body='Io', units = 'radec', wcs=w)
        
# Look up velocity of NH, for stellar aberration
        
        abcorr = 'LT+S'
        frame = 'J2000'
        st,ltime = cspice.spkezr('New Horizons', et, frame, abcorr, 'Sun') # Get velocity of NH 
        vel_sun_nh_j2k = st[3:6]
        
# Correct stellar RA/Dec for stellar aberration

        radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
        radec_stars_abcorr = hbt.correct_stellab(radec_stars, vel_sun_nh_j2k) # Store as radians

# Convert ring RA/Dec for stellar aberration

        radec_ring1        = np.transpose(np.array((ra_ring1,dec_ring1)))
        radec_ring1_abcorr = hbt.correct_stellab(radec_ring1, vel_sun_nh_j2k) # radians
        radec_ring2        = np.transpose(np.array((ra_ring2,dec_ring2)))
        radec_ring2_abcorr = hbt.correct_stellab(radec_ring2, vel_sun_nh_j2k) # radians
        
# Convert RA/Dec values back into pixels
        
        x_stars,        y_stars          = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)        
        x_stars_abcorr, y_stars_abcorr   = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
        x_ring1_abcorr, y_ring1_abcorr   = w.wcs_world2pix(radec_ring1_abcorr[:,0]*r2d, radec_ring1_abcorr[:,1]*r2d, 0)
        x_ring2_abcorr, y_ring2_abcorr   = w.wcs_world2pix(radec_ring2_abcorr[:,0]*r2d, radec_ring2_abcorr[:,1]*r2d, 0)

        points_stars        = np.transpose((x_stars, y_stars))
        points_stars_abcorr = np.transpose((x_stars_abcorr, y_stars_abcorr))

# Read the image file from disk

        image_polyfit = hbt.get_image_nh(t['Filename'][index_image], frac_clip = 1.,  
                                     bg_method = 'Polynomial', bg_argument = 4)
        image_raw     = hbt.get_image_nh(t['Filename'][index_image], frac_clip = 0.9, 
                                     bg_method = 'None')

# Use DAOphot to search the image for stars. It works really well.

        points_phot = find_stars(image_polyfit)
        
# Now look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
# I call this 'opnav'. It is returned in order (y,x) because that is what imreg_dft uses, even though it is a bit weird.
#
# For this, I can use either abcorr stars or normal stars -- whatever I am going to compute the offset from.        

        (dy_opnav, dx_opnav) = calc_offset_points(points_phot, points_stars, np.shape(image_raw), plot=False)

# Save the newly computed values to variables that we can access externally
# For the star locations, we can't put an array into an element of an astropy table.
# But we *can* put a string into astropy table! Do that: wrap with repr(), unwrap with eval('np.' + ).
       
        self.t_group['dx_opnav'][self.index_image] = dx_opnav
        self.t_group['dy_opnav'][self.index_image] = dy_opnav
        
        self.t_group['x_pos_star_cat'][self.index_image] = repr(x_stars_abcorr)
        self.t_group['y_pos_star_cat'][self.index_image] = repr(y_stars_abcorr)

#        self.t_group['x_pos_star_cat'][self.index_image] = repr(x_stars) # For test purposes, ignore the abcorr
#        self.t_group['y_pos_star_cat'][self.index_image] = repr(y_stars)
        
        self.t_group['x_pos_star_image'][self.index_image] = repr(points_phot[:,0])
        self.t_group['y_pos_star_image'][self.index_image] = repr(points_phot[:,1])
        self.t_group['is_navigated'][self.index_image] = True

        self.t_group['x_pos_ring1'][self.index_image] = repr(x_ring1_abcorr)
        self.t_group['y_pos_ring1'][self.index_image] = repr(y_ring1_abcorr)
        self.t_group['x_pos_ring2'][self.index_image] = repr(x_ring2_abcorr)
        self.t_group['y_pos_ring2'][self.index_image] = repr(y_ring2_abcorr)
        
        print "Opnav computed: {dx,dy}_opnav = " + repr(dx_opnav) + ', ' + repr(dy_opnav)

        print 'ra_stars      : ' + repr(ra_stars*r2d) + ' deg'
               
        # Now that we have navigated it, replot the image!
       
        self.plot_image()
        
        return 0

### THIS STUFF BELOW HERE IS NOT USED -- KEPT FOR HISTORICAL PURPOSES ONLY?? ###
        

######## TEST WCS OFFSETTING ##########

# Now convert the dx / dy offset (from opnav) into an RA/Dec offset.
# Put the new correct center position into WCS, so plots from now on will be correct.

        do_offset_wcs = False
        
        if (do_offset_wcs):
            m = w.wcs.piximg_matrix
            w.wcs.crval[0] -= (dy_opnav * (m[0,0] + m[1,0])) # RA. I am sort of guessing about these params, but looks good.
            w.wcs.crval[1] -= (dx_opnav * (m[1,1] + m[0,1])) # Dec
    
            ax1 = figs.add_subplot(1,2,2, projection=w) # nrows, ncols, plotnum. Returns an 'axis'
            fig1 = plt.imshow(np.log(image_polyfit))
            # Get the x and y pixel positions from WCS. The final argument ('0') refers to row vs. column order
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)       
            
            plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', 
                              markersize=12)
    
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)      
            x_stars_abcorr, y_stars_abcorr = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
            x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)
    
    
            plt.plot(x_stars, y_stars, marker='o', ls='None', color='lightgreen', mfc = 'None', mec = 'red', \
                              ms=12, mew=2, label='Cat')
            plt.plot(x_stars_abcorr, y_stars_abcorr, marker='o', ls='None', mfc = 'None', mec = 'green', \
                              ms=12, mew=2, label = 'Cat, abcorr')
            plt.plot(x_ring_abcorr, y_ring_abcorr, marker='o', ls='-', color='blue', mfc = 'None', mec = 'blue', 
                              ms=12, mew=2, 
                     label = 'Ring, LT+S')
                    
            plt.xlim((0,1000))
            plt.ylim((0,1000))
            plt.title('WCS, center = ' + repr(w.wcs.crval))

##########
# INIT CLASS
##########

    def __init__(self, master):

# Set some default values

        self.filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

        dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

        file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

# Set the sizes of the plots -- e.g., (7,7) = large square

        self.figsize_1 = (7,7) # Setting this larger and larger still leaves whitespace...
        self.figsize_2 = (7,3)
        self.figsize_3 = (7,4)
        
        option_bg_default = 'Grp Num Frac Pow'
        entry_bg_default  = '4' # Default polynomial order
        index_group_default = 2 # Jupiter ring phase curve
        index_image_default = 0 # Image number within the group

# Do some general code initialization

        hbt.set_plot_defaults()
        
        self.bgcolor = '#ECECEC'    # ECECEC is the 'default' background for a lot of the ttk widgets.
        
        (r2d, d2r) = (hbt.r2d, hbt.d2r)

# Start up SPICE

        cspice.furnsh(file_tm)

# Check if there is a pickle save file found. If it is there, go ahead and read it.

        if os.path.exists(self.filename_save):
            print "Loading file: " + self.filename_save
            self.load(verbose=False) # This will load self.t
            t = self.t

        else:    

# Find and process all of the FITS header info for all of the image files
# 't' is our main data table, and has all the info for all the files.
        
            t = hbt.get_fits_info_from_files_lorri(dir_images)

# Define new columns for the main table. 
# For several of these, we want to stick an array into the table cell. We do this by defining it as a string.
         
            t['dx_opnav']  = 0 # The measured offset between pos_star_cat and pos_star_image
            t['dy_opnav']  = 0
            t['dx_offset'] = 0 # Additional user-requested offset between the two
            t['dy_offset'] = 0
            t['bg_method'] = option_bg_default # None, Previous, Next, Polynomial, Median, "Grp Num Frac Pow", etc.
            t['bg_argument'] = entry_bg_default # A number (if 'Polynomial'). A range (if 'Median')
            t['Comment']  = 'Empty comment'  # Blank -- I'm not sure how to init its length if needed
            t['is_navigated'] = False  # Flag
            t['x_pos_star_cat']   = np.array(t['Format'],dtype='S5000')   # 1 x n array, pixels, with abcorr
            t['y_pos_star_cat']   = np.array(t['Format'],dtype='S5000')   # 1 x n array, pixels, with abcorr            
            t['x_pos_star_image'] = np.array(t['Format'],dtype='S5000') # 1 x n array, pixels, with abcorr
            t['y_pos_star_image'] = np.array(t['Format'],dtype='S5000') # 1 x n array, pixels, with abcorr
            t['x_pos_ring1']      = np.array(t['Format'],dtype='S5000')
            t['y_pos_ring1']      = np.array(t['Format'],dtype='S5000')
            t['x_pos_ring2']      = np.array(t['Format'],dtype='S5000')
            t['y_pos_ring2']      = np.array(t['Format'],dtype='S5000')
                      
            self.t = t
        
        # Get a list of unique observation descriptions (i.e., 'groups')
        
        self.groups = astropy.table.unique(t, keys=(['Desc']))['Desc']
                
        # Set initial location to first file in first group
        
        self.index_group = index_group_default
        self.index_image  = index_image_default
        
        # Extract the one group we are looking at 
        # NB: This creates an entirely new table -- *not* a view into the original table. 
        # If we modify values here, we need to explicitly write the changes back to the original.
        
        self.groupmask = t['Desc'] == self.groups[self.index_group]
        self.t_group = t[self.groupmask]
        
        self.num_images_group = np.size(self.t_group)

# Now set up a few variables for the current image only. These are not stored in the main
# table, because they are regenerated / reloaded every time, and they are potentially very large.

        self.planes             = 0         # Backplanes will be loaded into this variable
        self.has_backplane      = False     # Need to update this every time we load a new image file.
        self.profile_radius     = np.zeros((1))
        self.profile_azimuth    = np.zeros((1))
        self.bins_radius        = np.zeros((1))
        self.bins_azimuth       = np.zeros((1))
        
        self.image              = np.zeros((1)) # Save a copy of the current image (including subtraction, filtering, etc.)

        self.do_autoextract     = 1   # Flag. But make it 1/0, not True/False, as per ttk.        
        radii                   = cspice.bodvrd('JUPITER', 'RADII')
        self.rj                 = radii[0]  # Jupiter radius, polar, in km. Usually 71492.
        
# Now define and startup the widgets

# Create a container

#        frame = ttk.Frame(master, width=500)
                
        self.path_settings = "/Users/throop/python/salt_interact_settings/"

# Create the sliders, for navigation offset
        
        self.slider_offset_dx  = ttk.Scale(master, from_=-100, to=100, orient=Tkinter.HORIZONTAL, 
                                   command=self.set_offset) # Slider dx offset
        self.slider_offset_dy  = ttk.Scale(master, from_=-100, to=100, orient=Tkinter.VERTICAL, 
                                   command=self.set_offset) # Slider dy offset

# Define labels

        self.var_label_info = Tkinter.StringVar()
        self.var_label_info.set("init")
        self.label_info = ttk.Label(master, textvariable = self.var_label_info)
        
        self.var_label_status_io = Tkinter.StringVar()
        self.var_label_status_io.set('IO STATUS')
        self.label_status_io = ttk.Label(master, textvariable = self.var_label_status_io)       # For the 'save status' text
                   
# Create the buttons
   
        self.button_prev =     ttk.Button(master, text = 'Prev',     command=self.select_image_prev)
        self.button_next =     ttk.Button(master, text = 'Next',     command=self.select_image_next)
        self.button_plot   =   ttk.Button(master, text = 'Plot',     command=self.plot_image)
        self.button_extract =  ttk.Button(master, text = 'Extract',  command=self.extract_profiles)
        self.button_navigate = ttk.Button(master, text = 'Navigate', command=self.navigate)
        self.button_repoint =  ttk.Button(master, text = 'Repoint',  command=self.repoint)
        self.button_ds9 =      ttk.Button(master, text = 'Open in DS9', command=self.open_external)
        self.button_copy =     ttk.Button(master, text = 'Copy name to clipboard', command=self.copy_name_to_clipboard)
        self.button_quit =     ttk.Button(master, text = 'Quit',     command = self.quit)
        self.button_load =     ttk.Button(master, text = 'Load',     command = self.load)
        self.button_save =     ttk.Button(master, text = 'Save',     command = self.save_file)

        self.var_checkbutton_extract = Tkinter.IntVar()
        self.var_checkbutton_extract.set(self.do_autoextract)        # Values are in fact 1,0  not  True,False. Ugh. 
        self.checkbutton_extract = ttk.Checkbutton(master, text = 'Auto-Extract', command = self.handle_checkbutton, 
                                                   var = self.var_checkbutton_extract)

        
# Create controls for the background subtraction
# NB: Sometimes the current value for OptionMenu is left as blank. If this happens, need to restart spyder.
      
        self.var_option_bg = Tkinter.StringVar()
        self.var_option_bg.set("Polynomial Order:")
        
        self.label_bg =         ttk.Label(master, text='Background Method')
        self.option_bg =        ttk.OptionMenu(master, self.var_option_bg,\
                                    option_bg_default,  # First argument is the default
                                    "Polynomial", "Previous", "Next", "Median Range", "None", "Grp Num Frac Pow", \
                                     command = self.select_bg)
#        self.option_bg.config(width=25)
        self.entry_bg=          ttk.Entry(master, width=12)
        self.entry_bg.insert(0, entry_bg_default) # Set the value (e.g., order = "4")

                
# Create the Entry boxes (ie, single-line text)
# Note that for .insert(num, text), 'num' is in format line.character. Line starts at 1; char starts at 0.
# For both .text and .entry .
        
        self.entry_comment = ttk.Entry(master)
        self.entry_comment.insert(0, 'Comment')

# Create the text box to stuff the header into
# Since we are just displaying text, I'm not sure if this is the right widget to use
# http://www.tkdocs.com/tutorial/text.html
        
        self.text_header = Tkinter.Text(master, height=15)

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

        self.lbox_groups = Tkinter.Listbox(master, height=5, selectmode = "browse",
                                           width=25, exportselection=False)
        for item in self.groups.data:
            self.lbox_groups.insert("end", item)

        self.lbox_groups.bind('<<ListboxSelect>>', self.select_group) # Define the event handler for 
# Listbox is tied to an event handler, *not* to the normal command= thing.

# Create the listbox for files

        files_short = self.t_group['Shortname']
        num_files = np.size(files_short)

        self.lbox_files = Tkinter.Listbox(master, height=20, selectmode = "browse",
                                            width = 30, exportselection=False)

# Populate the listbox for files

        for i in range(np.size(files_short)):
            
             s = '{0:3}.   {1}   {2:.3f}   {3}  {4}'.format( \
                 repr(i), 
                 files_short[i], 
                 self.t_group['Exptime'][i], 
                 self.t_group['Format'][i],
                 self.t_group['UTC'][i]
                 
                 )
             self.lbox_files.insert("end", s)
                  
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
     
        self.fig1 = Figure(figsize = (7,7))
        junk = self.fig1.set_facecolor(self.bgcolor)

        self.ax1 = self.fig1.add_subplot(1,1,1, 
                                    xlabel = 'X', ylabel = 'Y', 
                                    label = 'Image',
                                    xlim = (0,1023),
                                    ylim = (1023,0)) # Return the axes
 
        self.canvas1 = FigureCanvasTkAgg(self.fig1,master=master)
        self.canvas1.show()
        self.ax1.imshow(hbt.get_image_nh('')) # Put up an empty frame, if file = ''
        
# Set up Plot 2 : Radial profile
        
        self.fig2 = Figure(figsize = (7,3))
        junk = self.fig2.set_facecolor(self.bgcolor)

        self.ax2 = self.fig2.add_subplot(1,1,1, 
                                    xlabel = 'Radius', ylabel = 'Intensity', 
                                    label = 'Plot') # Return the axes
 
        self.canvas2 = FigureCanvasTkAgg(self.fig2,master=master)
        self.canvas2.show()  

# Set up Plot 3 : Azimuthal profile
        
        self.fig3 = Figure(figsize = (7,4))
        junk = self.fig3.set_facecolor(self.bgcolor)

        self.ax3 = self.fig3.add_subplot(1,1,1, 
                                    xlabel = 'Azimuth', ylabel = 'Intensity', 
                                    label = 'Plot') # Return the axes
 
        self.canvas3 = FigureCanvasTkAgg(self.fig3,master=master)
        self.canvas3.show()  

# Set up Plot 4 : Debugging / Misc window
     
        self.fig4 = Figure(figsize = (4,4))
        junk = self.fig4.set_facecolor(self.bgcolor)

        self.ax4 = self.fig4.add_subplot(1,1,1, 
                                    xlabel = 'X', ylabel = 'Y', 
                                    label = 'Image',
                                    xlim = (0,1023),
                                    ylim = (1023,0)) # Return the axes
 
        self.canvas4 = FigureCanvasTkAgg(self.fig4,master=master)
        self.canvas4.show()
        plot4 = self.ax4.imshow(hbt.get_image_nh('')) # Put up an empty frame, if file = ''
        self.fig4.colorbar(plot4)
        
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

        self.button_navigate.        grid(row=8,  column=11)
        self.button_plot.            grid(row=8,  column=13)
        self.button_extract.         grid(row=8, column=14)
    
        self.button_prev.            grid(row=9,  column=11, columnspan=1, sticky='we')
        self.button_next.            grid(row=9,  column=13, columnspan=1, sticky='we')
        self.checkbutton_extract.    grid(row=9,  column=14)

        self.button_load.            grid(row=10,  column=11)
        self.button_save.            grid(row=10,  column=13)
        self.label_status_io.        grid(row=10,  column=14) 

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

#        self.save_settings()  # Small bug: for the first image, t[] is set, but not the gui, so plot is not made properly.
        
        self.load_backplane()
        
        self.refresh_gui()

##########
# Change Group
##########

    def select_group(self, event):

        print "select_group"

        index = (self.lbox_groups.curselection())[0]

        print "index = " + repr(index)
        
# First change the image # back to zero (or anything, but 0 is safest, and works even if currently at 0).
# This is necessary to save the state of the image (extraction params, etc.)

        self.index_image_new = 0
        self.change_image(refresh=False) # Save all parameters. 
        
# Now look up the new group name. Extract those elements.        
          
        name_group = self.groups[index]
        
        self.groupmask = self.t['Desc'] == name_group
        
        self.t_group = self.t[self.groupmask]

        self.num_images_group = np.size(self.t_group)

        print "Group #{} = {} of size {} hit!".format(index, name_group, self.num_images_group)

# Create the listbox for files, and load it up

        files_short = self.t_group['Shortname']

        self.lbox_files.delete(0, "end")
 
        for i in range(np.size(files_short)):
            
             s = '{0:3}.   {1}   {2:.3f}   {3}  {4}'.format( \
                 repr(i), 
                 files_short[i], 
                 self.t_group['Exptime'][i], 
                 self.t_group['Format'][i],                 
                 self.t_group['UTC'][i])

             self.lbox_files.insert("end", s)
                 
# Then set the group number, and refresh screen.
        
        self.index_group = index
        self.refresh_gui()

##########
# Save settings. Take all settings from the GUI, and put them into the proper table variables
##########
# This does *not* save to disk. It saves the Tk widget settings into variables, which can then be saved later.
 
    def save_settings(self):

#        print "save_settings"
        
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
# index_image has the current image (within the current group)
# index_image_new has the new image we want to show (within the current group)
# This does not change the group -- just the image      

    def change_image(self, refresh=True):
   
#        print "change_image"

        # Save current GUI settings into variables
        
        self.save_settings()
        
# Update the index, to show that new image is active
        
        self.index_image = self.index_image_new

# Load the new backplane, for the new image. Does not do any plotting.

        self.load_backplane()
        
# Load the new image and header, and show them
        
        self.refresh_gui()
        
# Extract radial / azimuthal profiles

        if (self.do_autoextract == 1):
            self.extract_profiles()        

##########
# Select new image, based on user click
##########
# This is called when user clicks on a new image number in list.

    def select_image(self, event):

        print "select_image"
        
        index = (self.lbox_files.curselection())[0]
        name = self.t_group['Shortname'][index]
        print "selected = " + repr(index) + ' ' + name
        self.index_image_new = index
        
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
# Update UI for new settings. Also update image if requested.
##########
# A general-purpose routine that is called when we have clicked on a new image, or otherwise changed
# the image number. This updates all of the GUI elements to reflect the new image, based on 
# the internal state which is already correct.
#
# This does *not* refresh the image itself, unless refresh_image=True is passed.

    def refresh_gui(self, refresh_image = True):
         
        if (refresh_image):
            self.plot_image()
        
# Load and display the proper header

        filename = self.t_group['Filename'][self.index_image]
          
        header = hbt.get_image_header(filename, single=True)
        self.text_header.delete(1.0, "end")
        h2 = re.sub('\s+\n', '\n', header)  # Remove extra whitespace -- wraps around and looks bad
        self.text_header.insert(1.0, h2)

        self.refresh_statusbar()

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

##########
# Plot image
##########

    def plot_image(self):
        """
        Load the specified file from disk, looks up the correct plotting parameters, and plots it.
        """


#        print "plot_image"
        
        # First save the current settings in case they've not beeen saved yet.
        # Then redraw the screen, based on new settings.
        
#        self.index_image_new = self.index_image
#        self.change_image()

        t = self.t_group[self.index_image]  # Grab this, read-only, since we use it a lot.
                                            # We can reference table['col'][n] or table[n]['col'] - either OK
         
        filename = self.t_group['Filename'][self.index_image]
        filename = str(filename)

        # Clear the image

        self.ax1.clear()
        
        # Set the color map
        
        plt.rc('image', cmap='Greys_r')
        
# Plot the image itself.
# This will always reload the image from disk.

        method = self.var_option_bg.get()

        print "OptionMenu: method = " + method
        
#        image = np.zeros((1023, 1023))
        
        frac, poly = 0, 0
        
        if (method == 'Previous'):
            file_prev = self.t_group['Filename'][self.index_image-1]
            print "file =      " + filename
            print "file_prev = " + file_prev
            image_bg = hbt.get_image_nh(file_prev, frac_clip = 1.0, bg_method = 'None')
            image_fg = hbt.get_image_nh(filename,  frac_clip = 1.0, bg_method = 'None')
            image = image_fg - image_bg

        if (method == 'Next'):
            file_next = self.t_group['Filename'][self.index_image+1]
            image_bg = hbt.get_image_nh(file_next, frac_clip = 1.0, bg_method = 'None')
            image_fg = hbt.get_image_nh(filename,  frac_clip = 1.0, bg_method = 'None')
            image = image_fg - image_bg
            
        if (method == 'Median'): # XXX not working yet
            file_prev = self.t_group['Filename'][self.index_image-1]
            image_bg = hbt.get_image_nh(file_prev, frac_clip = 1.0, bg_method = 'None')
            image_fg = hbt.get_image_nh(filename,  frac_clip = 1.0, bg_method = 'None')
            image = image_fg - image_bg

        if (method == 'Polynomial'):
            image = hbt.get_image_nh(filename,     frac_clip = 1.0, bg_method = 'Polynomial',
                                                   bg_argument = self.entry_bg.get())
                         
        if (method == 'Grp Num Frac Pow'):  # This is the most common method, I think
            vars = self.entry_bg.get().split(' ')
            
            if (np.size(vars) == 0): # If no args passed, just plot the image
                poly = 0
                frac = 0
                image = hbt.get_image_nh(filename,    frac_clip = 1.0, bg_method = 'None')

            if (np.size(vars) == 1): # One variable: interpret as exponent
                poly = float(vars[0])
                frac = 0
                image = hbt.get_image_nh(filename,    frac_clip = 1.0, bg_method = 'None')
                image = image - hbt.sfit(image, poly)
                
            if (np.size(vars) == 2): # Two variables: interpret as group num and file num
                (grp, num) = vars
                frac = 1
                poly = 0
                
            if (np.size(vars)) == 3: # Three variables: interpret as group num, file num, fraction
                (grp, num, frac) = vars
                poly = 0
                
            if (np.size(vars) == 4): # Four variables: Group num, File num, Fraction, Exponent
                (grp, num, frac, poly) = vars
                 
            if int(np.size(vars)) in [2,3,4]:
               
                grp = int(grp)
                num = int(num)
                frac = float(frac)
                poly = int(poly)
                
                print "group={}, num={}, frac={}".format(grp, num, frac)
#            print "Group = {}, num{}, Name = {}".format(name_group, num, name)

                name_group = self.groups[grp]
                groupmask = self.t['Desc'] == name_group
                group_tmp = self.t[groupmask]
                filename_bg = group_tmp['Filename'][num]
                                        
                image_fg = hbt.get_image_nh(filename,    frac_clip = 0.9, bg_method = 'None')
                image_bg = hbt.get_image_nh(filename_bg, frac_clip = 0.9, bg_method = 'None')
                
                image = image_fg - float(frac) * image_bg                
                image = image - hbt.sfit(image, poly)
                             
        if (method == 'None'):
            image = hbt.get_image_nh(filename, frac_clip = 1.0, bg_method = 'None')

        # Scale image by removeing stars, CRs, etc.

        image = hbt.remove_brightest(image, 0.97)   # Remove positive stars
        
        if (frac > 0):
            image = -hbt.remove_brightest(-image, 0.97) # Remove negative stars also
        
        # Save the image where we can use it later for extraction. 
        # XXX We might want to change this later: save a non-scaled, or non-subtracted version. 
        # But this is a good start.

        self.image = image
        
        # Render the main LORRI frame
        # *** In order to get things to plot in main plot window, 
        # use self.ax1.<command>, not plt.<command>
        # ax1 is an instance of Axes, and contains most other methods (legend, imshow, plot, etc)
                             
        self.ax1.imshow(image)

        # If we have them, draw the stars on top

        if (self.t_group['is_navigated'][self.index_image]):

            x_pos_ring1 = eval('np.' + t['x_pos_ring1']) # Convert from string (which can go in table) to array
            y_pos_ring1 = eval('np.' + t['y_pos_ring1'])
            x_pos_ring2 = eval('np.' + t['x_pos_ring2'])
            y_pos_ring2 = eval('np.' + t['y_pos_ring2'])           
            
            dx = t['dx_opnav'] + self.slider_offset_dx.get()
            dy = t['dy_opnav'] + self.slider_offset_dy.get()
            
            self.ax1.plot(eval('np.' + t['x_pos_star_cat']) + t['dx_opnav'], 
                     eval('np.' + t['y_pos_star_cat']) + t['dy_opnav'], 
                     marker='o', ls='None', 
                     color='lightgreen', ms=12, mew=1, label = 'Cat Stars, OpNav')
                     
            self.ax1.plot(eval('np.' + t['x_pos_star_cat']), 
                     eval('np.' + t['y_pos_star_cat']), 
                     marker='o', ls='None', 
                     color='lightgreen', label = 'Cat Stars, WCS')

            self.ax1.plot(eval('np.' + t['x_pos_star_image']), 
                     eval('np.' + t['y_pos_star_image']), 
                     marker='o', ls='None', 
                     color='pink', label = 'DAOfind Stars')               
            
            DO_PLOT_RING_INNER = False
            DO_PLOT_RING_OUTER = True
            
            if (DO_PLOT_RING_OUTER):
                self.ax1.plot(x_pos_ring2, y_pos_ring2, marker='o', ls = '--',
                              label = 'Ring, OpNav only')

                self.ax1.plot(x_pos_ring2 + dx, y_pos_ring2 + dy, marker='o', ls = '--',
                              label = 'Ring, OpNav+User')
                    
            if (DO_PLOT_RING_INNER):
                self.ax1.plot(x_pos_ring1, y_pos_ring1, marker='o', color='green', ls = '-', \
                    ms=8, label='Ring, LT')

                self.ax1.plot(x_pos_ring1 + dx, y_pos_ring1 + dy, \
                    marker='o', color='purple', ls = '-', ms=8, label='Ring, LT, Shifted')

            
            self.ax1.legend()
            
        # Disable the tickmarks from plotting

        self.ax1.get_xaxis().set_visible(False)
        self.ax1.get_yaxis().set_visible(False)

        # Set image size (so off-edge stars are clipped, rather than plot resizing)

        self.ax1.set_xlim([0,1023])  # This is an array and not a tuple. Beats me, like so many things with mpl.
        self.ax1.set_ylim([1023,0])
  
#        self.ax1.Axes(fig, [0,0,1,1])
        
        # Draw the figure on the Tk canvas
        
        self.fig1.tight_layout() # Remove all the extra whitespace -- nice!
        self.canvas1.draw()
        
        print "finished drawing canvas"
        
        return 0

##########
# Extract ring profiles (both radial and azimuthal), and plot them
##########

    def extract_profiles(self):
        
        if (self.t_group['is_navigated'][self.index_image] == False):
  
            print "Not navigated -- skipping profile"
            return 0 # XXX for debugging -- turn this off since slow
        
        dx_total =  -( self.t_group['dx_opnav'][self.index_image] +  int(self.slider_offset_dx.get()) )
        dy_total =  -( self.t_group['dy_opnav'][self.index_image] +  int(self.slider_offset_dy.get()) )
        
        image_roll = np.roll(np.roll(self.image, dx_total, axis=1), dy_total, axis=0)
        
        radius = self.planes['Radius_eq'] / self.rj   # Radius stored as km, but do the calc in R_J
        azimuth = self.planes['Longitude_eq'] * hbt.r2d # Azimuth is stored as radians, but do the calc in degrees
        
        print "radius[100,100] = " + repr(radius[100,100])

        num_bins_azimuth = 360 # Might make this a user parameter later
        num_bins_radius  = 360

        plt.rc('image', cmap='jet')               # Default color table for imshow

        self.ax4.imshow( ( np.array(radius > 1.7) & np.array(radius < 1.8)) + 1* image_roll / np.max(image_roll))
        self.canvas4.show()
        
        # Extract and plot radial profile

        self.bins_radius = hbt.frange(1.6, 1.9, num_bins_radius)
        self.profile_radius = np.zeros(num_bins_radius)
        
        for i in range(num_bins_radius-1):
            
            is_good = np.array(radius > self.bins_radius[i])                  & np.array(radius < self.bins_radius[i+1])
                      
            self.profile_radius[i] = np.mean(image_roll[is_good])
        
        self.ax2.plot(self.bins_radius, self.profile_radius)
        fs = 20
#        self.ax2.set_xlim([0,100])  # This is an array and not a tuple. Beats me, like so many things with mpl.

        self.canvas2.show()

        # Extract and plot azimuthal profile
        
        # Create the azimuthal bins. 
        self.bins_azimuth = hbt.frange(-180,180,num_bins_azimuth+1)
        self.profile_azimuth = np.zeros(num_bins_azimuth+1)
                
        print
        for i in range(num_bins_azimuth-1):
            
            is_good = np.array(radius > 1.5)                  & np.array(radius < 2.2) & \
                      np.array(azimuth > self.bins_azimuth[i]) & np.array(azimuth < self.bins_azimuth[i+1])
                      
            self.profile_azimuth[i] = np.mean(image_roll[is_good])
#            print "{:<3}. {:<7} {:<7}".format(i, self.bins_azimuth[i], self.profile_azimuth[i])

        # Now we do some crazy logic to unwrap the azimuth, so that start will be in -180 .. 180, and end after that.
        # First, copy the azimuth and intensity, so we have two full loops in the array (720 deg)

        profile_azimuth_2      = np.concatenate((self.profile_azimuth[:-2], self.profile_azimuth[:-2]))
        bins_azimuth_2         = np.concatenate((self.bins_azimuth[:-2], self.bins_azimuth[:-2]+360))     
        
        # Now, search this for the first pattern of [nan, non-nan]. That will be the start of the valid data.
        
        i = np.array(range(np.size(profile_azimuth_2)-1)) # Just an index. Chop off last entry so we can use i+1
        
        is_start = np.isnan(profile_azimuth_2[i]) & np.logical_not(np.isnan(profile_azimuth_2[i+1]))
        bin_start = i[is_start][0]+1 # get first match
        
        # And then look for the end pattern: [non-nan, nan], starting *after* bin_start
        
        is_end = np.logical_not(np.isnan(profile_azimuth_2[i])) & np.isnan(profile_azimuth_2[i+1]) & np.array(i > bin_start)
        bin_end = i[is_end][0]

        # Now we have the proper indices for the start and end.
 
        az_start = bins_azimuth_2[bin_start]
        az_end   = bins_azimuth_2[bin_end]
        
        print "az = " + repr(bins_azimuth_2[bin_start]) + ' .. ' + repr(bins_azimuth_2[bin_end])
        
        self.ax3.plot(bins_azimuth_2, profile_azimuth_2)
        fs = 20

#        self.ax3.set_xlim([az_start-10, az_end+10])  # This is an array and not a tuple. Beats me, like so many things with mpl.
#        self.ax3.set_ylim([-5, 5])           # Hard-code this.... not sure what is best.
        self.canvas3.show()
        
        plot4 = self.ax4.imshow(azimuth)
#        plt.title('Lon_eq [deg]', fontsize=fs)
#        plt.xlim((0,1000))
#        plt.ylim((0,1000))
    
#        self.fig4.colorbar(plot4)  # Need to get this to updated itself. Not sure how.
        self.canvas4.show()
            
        return 0

##########
# Load backplane
##########
# Returns True if loaded successfully; False if not

    def load_backplane(self):

        dir_backplanes = '/Users/throop/data/NH_Jring/out/'
        file_backplane = dir_backplanes + self.t_group['Shortname'][self.index_image].replace('.fit', '_planes.pkl')

        if (os.path.isfile(file_backplane)):
            print 'File_backplane = ' + file_backplane 
            lun = open(file_backplane, 'rb')
            planes = pickle.load(lun)
            self.planes = planes # This will load self.t
            lun.close()
            return True

        else:
            print "Can't load backplane " + file_backplane
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
        
        print "bg method = " + self.var_option_bg.get()

        # Grab the GUI settings for the background, and save them into proper variables.

        self.t_group['bg_method'][self.index_group] = self.var_option_bg.get()
        self.t_group['bg_argument'][self.index_group] = self.entry_bg.get()

##########
# Copy current filename (short) to clipboard
##########

    def copy_name_to_clipboard(self):

        str = self.t_group['Shortname'][self.index_image]
        hbt.write_to_clipboard(str)
        print "{} copied to clipboard".format(str)
    
        
        return 0
        
##########
# Close app and quit
##########
# Define this with one argument if we press a button, or two if we trigger an event like a keypress
# Get error if called with more than it can handle. So we define it here with two.

    def quit_e(self, event):
        self.quit()        
    
    def quit(self):
        root.destroy() # Q: How does root. get automatically imported here? Not sure.
        root.quit()

    def handle_checkbutton(self):
        
        print "Button = " + repr(self.var_checkbutton_extract.get())
        self.do_autoextract = self.var_checkbutton_extract.get() # Returns 1 or 0
        
##########
# Load settings
##########
# Load them from the pre-set filename

    def load(self, verbose=True):
        
        lun = open(self.filename_save, 'rb')
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
        self.save_settings()
    
        # Write one variable to a file    
    
        lun = open(self.filename_save, 'wb')
        t = self.t
        pickle.dump(t, lun)
        lun.close()
        
        # Put up a message for the user. A dialog is too intrusive.
                
        self.var_label_status_io.set('SAVING...')
        print "Set saving..."
        time.sleep(0.5)
        self.var_label_status_io.set('')

#            tkMessageBox.showinfo("Saved", "File " + self.filename_save + " saved.")
        
#        tk
##########
# (p)revious image
##########
      
    def select_image_prev(self):
        
        self.index_image_new = (self.index_image - 1) % self.num_images_group
        self.change_image()

##########
# (n)ext image
##########
    
    def select_image_next(self):
          
        self.index_image_new = (self.index_image + 1) % self.num_images_group
        self.change_image()

##########
# Process slider positions
##########

    def set_offset(self, event): # Because this is an event, it it passed two arguments
        
# Get the slider positions, for the dx and dy nav offset positions, and put them into a variable we can use

        self.offset_dx = self.slider_offset_dx.get() # 
        self.offset_dy = self.slider_offset_dy.get() # 
        
        self.plot_image()
        
        self.refresh_statusbar()
        
        return 0
        
###########
# Now start the main app
###########

# Start up the widget

root = Tkinter.Tk()
app  = App(root)

# set the dimensions and position of the window


root.geometry('%dx%d+%d+%d' % (1750, 900, 2, 2))
root.configure(background='#ECECEC')                # ECECEC is the 'default' background for a lot of the ttk widgets, which
                                                    # I figured out using photoshop. Not sure how to change that.

#  window.attributes('-topmost', 1)
#  window.attributes('-topmost', 0)

os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

    
root.mainloop()  # This will call the __init__ function


    



