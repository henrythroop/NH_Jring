#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:54:05 2016

@author: throop
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

@author: throop

"""

# -*- coding: utf-8 -*-
"""
# Program display a GUI of NJ J-ring data to navigate and extract it.

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
import wcsaxes
import time
from scipy.interpolate import griddata

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
import warnings

# HBT imports

import hbt

# First we define any general-purpose functions, which are not part of the class/module.
# We can move these to a different file at some point.

#########
# Navigate the image
#########

dir_data = '/Users/throop/Dropbox/Data/NH_Jring/data/jupiter/level2/lor/all/'
file     = dir_data + 'lor_0034765323_0x630_sci_1.fit'
hdulist  = fits.open(file)
arr      = hdulist['PRIMARY'].data
header   = hdulist['PRIMARY'].header

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")    
    w = WCS(file)                # Look up the WCS coordinates for this frame
                                 # Otherwise it gives "FITSFixedWarning: 'unitfix': 'Changed units: 'DEG' -> 'deg'"
    
et = header['SPCSCET']
utc = 
           

center  = w.wcs.crval  # degrees. # crval is a two-element array of [RA, Dec], in degrees
DO_GSC1     = False    # Stopped working 2-Oct-2016
DO_GSC2     = True
DO_USNOA2   = False
        

# Stretch the image
stretch_percent = 90
stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales array to 5th .. 95th %ile. 

# Display it for fun

plt.imshow(stretch(arr))

# Load
def navigate(self):        

    t = self.t_group # We use this a lot, so make it shorter
    d2r = hbt.d2r
    r2d = hbt.r2d
    
    index_image = self.index_image       
    
# Now look up positions of stars in this field, from a star catalog

    w = WCS(t['Filename'][index_image])                  # Look up the WCS coordinates for this frame
    et = t['ET'][index_image]

    print('ET[i] =  ' + repr(et))
    print('UTC[i] = ' + repr(t['UTC'][index_image]))
    print('crval[i] = ' + repr(w.wcs.crval))              # crval is a two-element array of [RA, Dec], in degrees
    
    center  = w.wcs.crval  # degrees
    DO_GSC1     = False    # Stopped working 2-Oct-2016
    DO_GSC2     = True
    DO_USNOA2   = False
    
    if (DO_GSC1):
        name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
        radius_search = 0.15
        stars = conesearch.conesearch(w.wcs.crval, radius_search, cache=False, catalog_db = name_cat)
        ra_stars  = np.array(stars.array['RAJ2000'])*d2r # Convert to radians
        dec_stars = np.array(stars.array['DEJ2000'])*d2r # Convert to radians
#            table_stars = Table(stars.array.data)

    if (DO_GSC2):
        name_cat = u'Guide Star Catalog v2 1'
#            name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating

#            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
        from astropy.utils import data
        
        with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                                       # The default is 3 seconds, and that times out often.
            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)

        ra_stars  = np.array(stars.array['ra'])*d2r # Convert to radians
        dec_stars = np.array(stars.array['dec'])*d2r # Convert to radians

        mag       = np.array(stars.array['Mag'])
        
        print("Stars downloaded: {}; mag = {} .. {}".format(np.size(mag), np.nanmin(mag), np.nanmax(mag)))
        print("RA = {} .. {}".format(np.nanmin(ra_stars)*r2d, np.nanmax(ra_stars)*r2d))
        
        # Now sort by magnitude, and keep the 100 brightest
        # This is because this GSC catalog is huge -- typically 2000 stars in LORRI FOV.
        # We need to reduce its size to fit in our fixed astropy table string length.

        num_stars_max = 100            
        order = np.argsort(mag)
        order = np.array(order)[0:num_stars_max]

        ra_stars = ra_stars[order]
        dec_stars = dec_stars[order]
  
#            table_stars = Table(stars.array.data)

    if (DO_USNOA2):        
        name_cat = u'The USNO-A2.0 Catalogue (Monet+ 1998) 1' # Works but gives stars down to v=17; I want to v=13 
        stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
        table_stars = Table(stars.array.data)
        mask = table_stars['Bmag'] < 13
        table_stars_m = table_stars[mask]            

        ra_stars  = table_stars_m['RAJ2000']*d2r # Convert to radians
        dec_stars = table_stars_m['DEJ2000']*d2r # Convert to radians
    
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
    st,ltime = sp.spkezr('New Horizons', et, frame, abcorr, 'Sun') # Get velocity of NH 
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

    image_polyfit = hbt.read_lorri(t['Filename'][index_image], frac_clip = 1.,  
                                 bg_method = 'Polynomial', bg_argument = 4)
    image_raw     = hbt.read_lorri(t['Filename'][index_image], frac_clip = 0.9, 
                                 bg_method = 'None')

# Use DAOphot to search the image for stars. It works really well.

    points_phot = hbt.find_stars(image_polyfit)
    
# Now look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
# I call this 'opnav'. It is returned in order (y,x) because that is what imreg_dft uses, even though it is a bit weird.
#
# For this, I can use either abcorr stars or normal stars -- whatever I am going to compute the offset from.        

    (dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot, points_stars, np.shape(image_raw), plot=False)

# Save the newly computed values to variables that we can access externally
# For the star locations, we can't put an array into an element of an astropy table.
# But we *can* put a string into astropy table! Do that: wrap with repr(), unwrap with eval('np.' + ).
   
    self.t_group['dx_opnav'][self.index_image] = dx_opnav
    self.t_group['dy_opnav'][self.index_image] = dy_opnav
    
    self.t_group['x_pos_star_cat'][self.index_image] = hbt.reprfix(x_stars_abcorr)
    self.t_group['y_pos_star_cat'][self.index_image] = hbt.reprfix(y_stars_abcorr)
    
#        self.t_group['x_pos_star_cat'][self.index_image] = repr(x_stars) # For test purposes, ignore the abcorr
#        self.t_group['y_pos_star_cat'][self.index_image] = repr(y_stars)
            
    self.t_group['x_pos_star_image'][self.index_image] = hbt.reprfix(points_phot[:,0]).replace(' ', '') # Shorten it
    self.t_group['y_pos_star_image'][self.index_image] = hbt.reprfix(points_phot[:,1]).replace(' ', '')
    self.t_group['is_navigated'][self.index_image] = True             # Set the flag saying we have navigated image

    self.t_group['x_pos_ring1'][self.index_image] = hbt.reprfix(x_ring1_abcorr)
    self.t_group['y_pos_ring1'][self.index_image] = hbt.reprfix(y_ring1_abcorr)
    self.t_group['x_pos_ring2'][self.index_image] = hbt.reprfix(x_ring2_abcorr)
    self.t_group['y_pos_ring2'][self.index_image] = hbt.reprfix(y_ring2_abcorr)
    
#        self.t_group['x_pos_bodies'][self.index_image] = repr(x_pos_bodies)
#        self.t_group['y_pos_bodies'][self.index_image] = repr(y_pos_bodies)
#        self.t_group['ra_bodies'][self.index_image]  = repr(ra_bodies)
#        self.t_group['dec_bodies'][self.index_image] = repr(dec_bodies)
#        self.t_group['name_bodies'][self.index_image] = repr(name_bodies)

#        print("Opnav computed: {dx,dy}_opnav = " + repr(dx_opnav) + ', ' + repr(dy_opnav))
#
#        print('ra_stars      : ' + repr(ra_stars*r2d) + ' deg')
           
    # Now that we have navigated it, replot the image!
           
    return 0

    """
    Load current image from disk
    """
    
    def load_image(self):

#        print("load_image()")

# autozoom: if set, and we are loading a 4x4 image, then scale it up to a full 1x1.

        file = self.t_group[self.index_image]['Filename']                    
        self.image_raw = hbt.read_lorri(file, frac_clip = 1., bg_method = 'None', autozoom=True)
        print("Loaded image: " + file)

# Load the backplane as well. We need it to flag satellite locations during the extraction.
        
        DO_LOAD_BACKPLANE = True
        if DO_LOAD_BACKPLANE:
            self.load_backplane()

#        print("Loaded backplane.")
                                                


#==============================================================================
# Apply image processing to image. Stray light, polynomial subtraction, etc.
#==============================================================================

    def process_image(self):
        
#        print("process_image()")
        
        method = self.var_option_bg.get()
        argument = self.entry_bg.get()
        
        self.image_processed = hbt.nh_jring_process_image(self.image_raw, \
                                                          method, argument, self.index_group, self.index_image)

        # Remove cosmic rays
        
        self.image_processed = hbt.decosmic(self.image_processed)
        
#==============================================================================
# Plot image
#==============================================================================

    def plot_image(self):
        """
        Plot the image. It has already been loaded and processed.
        This is an internal routine, which does not process.	
        This plots the image itself, *and* the objects, if navigated.							
        """
                
        # Clear the image

        self.ax1.clear()
        
        # Set the color map
        
        plt.rc('image', cmap='Greys_r')
        
# Plot the image itself.
        
        # Render the main LORRI frame
        # *** In order to get things to plot in main plot window, 
        # use self.ax1.<command>, not plt.<command>
        # ax1 is an instance of Axes, and I can call it with most other methods (legend(), imshow(), plot(), etc)

        stretch = astropy.visualization.PercentileInterval(self.stretch_percent)  # PI(90) scales array to 5th .. 95th %ile. 

# Get the satellite position mask, and roll it into position
                       
        self.ax1.imshow(stretch(self.image_processed))
        
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
                
        if (self.t_group['is_navigated'][self.index_image]):
            self.plot_objects()									

        
        plt.imshow(stretch(self.image_processed))
        plt.show()
        
        return 0

##########
# Plot Objects
##########

    def plot_objects(self):
        """
           Plot rings and stars, etc, if we have them.
      """
        
        # If we have them, draw the stars on top

        t = self.t_group[self.index_image]  # Grab this, read-only, since we use it a lot.
                                            # We can reference table['col'][n] or table[n]['col'] - either OK
        w = WCS(t['Filename'])                  # Look up the WCS coordinates for this frame
         
        filename = self.t_group['Filename'][self.index_image]
        filename = str(filename)
        
#        print("is_navigated: " + repr(self.t_group['is_navigated'][self.index_image])  )                      
								
        if (self.t_group['is_navigated'][self.index_image]):

            # Remove all of the current 'lines' (aka points) from the plot. This leaves the axis and the image
            # preserved, but just removes the lines. Awesome. Using ax1.cla() will clear the entire 'axis', including image.
            # Q: Does this remove legend? Not sure.
		
            lines = self.ax1.get_lines()
            for line in lines:
                line.remove()		# Not vectorized -- have to do it one-by-one

# Get ring position
												
            x_pos_ring1 = eval('np.' + t['x_pos_ring1']) # Convert from string (which can go in table) to array
            y_pos_ring1 = eval('np.' + t['y_pos_ring1'])
            x_pos_ring2 = eval('np.' + t['x_pos_ring2'])
            y_pos_ring2 = eval('np.' + t['y_pos_ring2'])           

# Get the user offset position
            
            dx = t['dx_opnav'] + self.slider_offset_dx.get()
            dy = t['dy_opnav'] + self.slider_offset_dy.get()

# Plot the stars -- catalog, and DAO

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

            # Get position of satellites

            name_bodies = np.array(['Metis', 'Adrastea', 'Thebe', 'Amalthea', 'Io'])        
            x_bodies,  y_bodies   = hbt.get_pos_bodies(t['ET'], name_bodies, units='pixels', wcs=w)
            ra_bodies, dec_bodies = hbt.get_pos_bodies(t['ET'], name_bodies, units='radec', wcs=w)

            # Plot satellites
              
            self.ax1.plot(x_bodies+dx, y_bodies+dy, marker = '+', color='red', markersize=20, linestyle='none')
            
            # Plot the ring
                
            DO_PLOT_RING_INNER = False
            DO_PLOT_RING_OUTER = True
            
            if (DO_PLOT_RING_OUTER):
                self.ax1.plot(x_pos_ring2, y_pos_ring2, marker='o', color = 'blue', ls = '--',
                              label = 'Ring, OpNav only')

                self.ax1.plot(x_pos_ring2 + dx, y_pos_ring2 + dy, marker='o', color = 'lightblue', ls = '--',
                              label = 'Ring, OpNav+User')
                    
            if (DO_PLOT_RING_INNER):
                self.ax1.plot(x_pos_ring1, y_pos_ring1, marker='o', color='green', ls = '-', \
                    ms=8, label='Ring, LT')

                self.ax1.plot(x_pos_ring1 + dx, y_pos_ring1 + dy, \
                    marker='o', color='purple', ls = '-', ms=8, label='Ring, LT, Shifted')

            self.legend = self.ax1.legend()  # Draw legend. Might be irrel since remove() might keep it; not sure.

            self.canvas1.draw()

        




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
    
        f_solar_1au_si     = 1.77                 # W/m2/nm. At 600 nm. 
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
        


    



