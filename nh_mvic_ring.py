# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 22:55:34 2016

@author: throop
"""

dir = '/Users/throop/Data/NH_MVIC_Ring/'
file = 'mvic_d305_sum_mos_v1.fits'

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
import astropy.visualization
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools   import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy     import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils   import daofind
import wcsaxes
import time
from scipy.interpolate import griddata

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

###

import hbt

dir     = '/Users/throop/Data/NH_MVIC_Ring' 
dir_out = dir + '/out'

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

cspice.furnsh(file_tm) # Start up SPICE

#==============================================================================
# Initialize constants
#==============================================================================

# '-new-image' is added by astrometry.net
# '_fixed'     is added by me to indicate that I have fixed the header (that it, added SCET field)
#files = ['mvic_d305_sum_mos_v1.fits', 'ringdep_mos_v1.fits']

#file = dir + '/mvic_d305_sum_mos_v1.fits'

sequence        = 'D305'
#sequence        = 'A_RINGDEP_01'  # Departure imaging, closest MVIC image of the whole system
#sequence        = 'D202'
 
DO_FIX_FITS     = False
DO_ANALYZE      = True

nbins_radius = 100
#nbins_radius = 1000

#PIXSIZE =              13.0000 /Pixel size in microns                           
#READNOI =              30.0000 /Readnoise in Electrons                          
#GAIN    =              58.6000 /Gain in Electrons/DN                            
#PIXFOV  =              19.8065 /Plate scale in microrad/pix    

#PSOLAR  = '2.5541E+14'         /(DN/s)/(erg/cm^2/s/Ang), Solar spectrum. Use for unresolved
#RSOLAR  = '100190.64'          /(DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum. Use for resolved sources.

exptime      = 10        # All MVIC Framing ring mosaics are 10 sec/exposure
rsolar       = 100190.64 # (DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum. Use for resolved sources.
pixfov       = 19.8065   # Plate scale in microrad/pix    
pix_sr       = (pixfov*(1e-6))**2



      = 1700 / (math.pi * 40.**2) * 1e7  # erg/cm2/sec at Pluto, roughly

#==============================================================================
# Initialize parameters for each of the possible MVIC mosaics we can analyze
#==============================================================================

# *_wcs.fits file are created by astrometry.net. Their default name is new-image.fits .
# *_header.fits indicates that I have added ET and any other necessary header info into the file. 
# *_pl.fits indicates that I have added backplanes.

if (sequence == 'D305'):
  file_wcs    = dir + '/mvic_d305_sum_mos_v1_wcs.fits' # Load the navigated image, with WCS  
  utc = '2015::305 00:00:00'  # Set the time 
  stretch = astropy.visualization.PercentileInterval(99)
  
if (sequence == 'A_RINGDEP_01'):
  file_wcs = dir + '/ringdep_mos_v1_wcs.fits'
  utc = '2015 Jul 15 18:50:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)

if (sequence == 'D202'):
  file_wcs = dir + '/mvic_d202_mos_v1_wcs.fits'
  utc = '2015::202 00:00:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
  
file_header    = file_wcs.replace('.fits', '_header.fits')       # File with fixed FITS ET info
file_header_pl = file_header.replace('.fits', '_pl.fits')     # File for backplanes

#==============================================================================
# Repair / edit the FITS files if needed
#==============================================================================
    
if (DO_FIX_FITS):

# Add a missing header field and write out
    
    et = cspice.utc2et(utc)
    hdulist = fits.open(file_wcs)
    im = hdulist['PRIMARY'].data    
    hdulist['PRIMARY'].header['SPCSCET'] =     (et, "[s past J2000] Spacecraft mid-obs time, TDB")  
    
    hdulist.writeto(file_header, clobber=True)
    print "Wrote new FITS file with fixed header: " + file_header

# Create the backplanes and write out

    planes = hbt.create_backplane(file_header, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')
    
#    lun = open(file_planes, 'wb')
#    pickle.dump(planes, lun)
#    lun.close()
    
#    hdulist[']
    print "Wrote backplane file: " + file_planes    

# Create backplaned FITS file.

#    hdu1 = fits.PrimaryHDU()
#    hdu2 = fits.ImageHDU()
#    new_hdul = fits.HDUList([hdu1, hdu2])
#    new_hdul.writeto('test.fits', clobber=True)

    type_out = 'float32'  # Write backplane in single precision, to save space.
    hdu_out = fits.HDUList()
    hdu_out.append(hdulist['PRIMARY'])
    hdu_out.append(fits.ImageHDU(planes['RA'].astype(type_out), name = 'RA'))
    hdu_out.append(fits.ImageHDU(planes['Dec'].astype(type_out), name = 'Dec'))
    hdu_out.append(fits.ImageHDU(planes['Phase'].astype(type_out), name = 'Phase'))
    hdu_out.append(fits.ImageHDU(planes['Longitude_eq'].astype(type_out), name = 'Longitude_eq'))
    hdu_out.append(fits.ImageHDU(planes['Radius_eq'].astype(type_out), name = 'Radius_eq'))
    
    hdu_out.writeto(file_header_pl, clobber=True)
    print "Wrote file: " + file_header_pl   

#==============================================================================
# Load the image + backplane
#==============================================================================
    
file = file_header_pl
hbt.figsize((6,28))

hdulist = fits.open(file)
im = hdulist['PRIMARY'].data
    
#lun = open(file_planes, 'rb')
#planes = pickle.load(lun)
#lun.close()

#radius = planes['Radius_eq']

radius    = hdulist['Radius_eq'].data
longitude = hdulist['Longitude_eq'].data

et = hdulist['PRIMARY'].header['SPCSCET'] # Get the ET from the file
utc = cspice.et2utc(et, 'C', 0)

w = WCS(file)

#==============================================================================
# Make a plot of the image + backplanes
#==============================================================================

plt.set_cmap('Greys_r')

hbt.figsize((20,40))
plt.subplot(1,3,1)
plt.imshow(stretch(im))
plt.title(sequence)
plt.gca().get_xaxis().set_visible(False)

plt.subplot(1,3,2)
plt.imshow(radius, cmap='plasma')
plt.title('Radius')
plt.gca().get_xaxis().set_visible(False)

plt.subplot(1,3,3)
plt.imshow(longitude, cmap='plasma')
plt.title('Longit')
plt.gca().get_xaxis().set_visible(False)

plt.show()

#==============================================================================
# Make a plot of radial profile
#==============================================================================
    
#hbt.figsize((5,4))
#plt.plot(bins_radius / 1000, flux_arr)
#plt.xlabel('Distance from Pluto [1000 km]')
#if (sequence == 'A_RINGDEP_01'):
#    plt.ylim((-0.1,0.1))
#plt.show()
    
#==============================================================================
# Make a labled plot with Pluto, Charon, Nix, Hydra, etc.
#==============================================================================

hbt.figsize((5,12))
plt.imshow(stretch(im))

offset_x, offset_y = 100,100  # Offset PCNHSK labels by this many pixels

pos_body_pix = {}  # Create a new dictionary to save positions of each body

name_body = ['Pluto', 'Charon', 'Nix', 'Hydra', 'Styx', 'Kerberos']
for name_body_i in name_body:
    vec,lt = cspice.spkezr(name_body_i, et, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    (junk,ra,dec) = cspice.recrad(vec_sc_targ) # Get the RA / Dec of the object
    x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
    plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
             markeredgecolor='red', markeredgewidth=2)
    plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='red', fontsize=15)

    pos_body_pix[name_body_i] = np.array([y_pix, x_pix])  # Save position of each body as a dictionary entry tuple
        
plt.title(sequence)
plt.xlim((0,np.shape(im)[1]))
plt.ylim((0,np.shape(im)[0]))
plt.show()

#plt.imshow(dist_pluto_pix)
#plt.xlim((0,np.shape(im)[1]))
#plt.ylim((0,np.shape(im)[0]))
#plt.title('Distance from Pluto')
#plt.show()
#
#plt.imshow(dist_charon_pix)
#plt.xlim((0,np.shape(im)[1]))
#plt.ylim((0,np.shape(im)[0]))
#plt.title('Distance from Charon')
#plt.show()

#==============================================================================
# Create new backplane set, with pixel distances from each body
#==============================================================================

dist_body_pix = {}

for name_body_i in name_body:
    
    (dist_body_pix_y,dist_body_pix_x) = np.meshgrid(\
                     np.array(range(np.shape(im)[1])) - int(pos_body_pix[name_body_i][1]),
                     np.array(range(np.shape(im)[0])) - int(pos_body_pix[name_body_i][0]))


    dist_body_pix[name_body_i] = np.sqrt(dist_body_pix_x**2 + dist_body_pix_y**2)

#==============================================================================
# Filter stars, CRs, and other outliers out of the image.
#==============================================================================

# Idea: look for pixels that are more that 2sigma away from mean of the image

std = np.std(im)
im_clean = im.copy()
is_outlier = im > (np.median(im) + 3*std)
im_clean[is_outlier] = np.median(im)

# Do it again

std = np.std(im_clean)
im_clean2 = im_clean.copy()
is_outlier = im_clean > (np.median(im_clean) + 3*std)
im_clean2[is_outlier] = np.median(im_clean)

#==============================================================================
# Measure the radial profile, from Pluto
#==============================================================================


bins_radius = hbt.frange(np.amin(radius), np.amax(radius), nbins_radius) # Set up radial bins, in km

dradius = bins_radius[1] - bins_radius[0]

bin_number_2d = np.copy(im)*0

flux_mean_arr   = np.zeros(nbins_radius)
flux_median_arr = np.zeros(nbins_radius)
flux_std_arr    = np.zeros(nbins_radius)
npix_arr        = np.zeros(nbins_radius)

flux_mean_clean_arr   = np.zeros(nbins_radius)
flux_median_clean_arr = np.zeros(nbins_radius)

for i in range(nbins_radius-1):
    is_good = np.array(radius > bins_radius[i]) \
            & np.array(radius < bins_radius[i+1]) \
            & (dist_body_pix['Pluto']  > 100) \
            & (dist_body_pix['Charon'] > 100) \
            & (im != 0)  # Mask out the zero pixels (e.g., between MVIC mosaic frames)
                         # Also need to mask out area near Pluto and near Charon
            
    bin_number_2d[is_good] = i          # For each pixel in the array, assign a bin number

    flux_mean_arr[i]         = np.mean(im[is_good])
    flux_mean_clean_arr[i]   = np.mean(im_clean2[is_good])
    flux_median_arr[i]       = np.median(im[is_good])
    flux_median_clean_arr[i] = np.median(im_clean2[is_good])
    flux_std_arr[i]          = np.std(im[is_good]) # stdev
    
    npix_arr[i] = np.sum(is_good) # Count of number of pixels in this radial bin


#==============================================================================
# Calc orbital distances (from Pluto) for each body
#==============================================================================

d_pluto_body = {}
mask_orbit = {}

if (sequence == 'D305'):
    d_radius = 4000  # Halfwidth to plot, of orbit, in km
  
if (sequence == 'A_RINGDEP_01'):
    d_radius = 100
  
if (sequence == 'D202'):
    d_radius = 500
    
# Look up the distance using SPICE

for name_body_i in name_body:
  (vec_body,junk) = cspice.spkezr(name_body_i, et, 'IAU_PLUTO', 'LT', 'Pluto')
  d_pluto_body[name_body_i] = cspice.vnorm(vec_body[0:3])
   
# Generate a pixel mask showing the orbit of each body

  mask_orbit[name_body_i] = \
    np.array(radius > (d_pluto_body[name_body_i] - d_radius)) & \
    np.array(radius < (d_pluto_body[name_body_i] + d_radius))

    
r_h = d_pluto_body['Hydra']  # Hydra orbital radius

mask_orbit['Hydra x 4'] = \
    np.array(radius > (r_h*4 - d_radius)) & np.array(radius < (r_h*4 + d_radius))
    
mask_orbit['Hydra x 10'] = \
    np.array(radius > (r_h*10 - d_radius)) & np.array(radius < (r_h*10 + d_radius))

mask_orbit['Hydra x 20'] = \
    np.array(radius > (r_h*20 - d_radius)) & np.array(radius < (r_h*20 + d_radius))

mask_orbit['Hydra x 40'] = \
    np.array(radius > (r_h*40 - d_radius)) & np.array(radius < (r_h*40 + d_radius))

mask_orbit['Hydra x 60'] = \
    np.array(radius > (r_h*60 - d_radius)) & np.array(radius < (r_h*60 + d_radius))
    
mask_orbit['Hydra x 80'] = \
    np.array(radius > (r_h*80 - d_radius)) & np.array(radius < (r_h*80 + d_radius))    
    
#==============================================================================
# Make plots of radial profile
#==============================================================================

hbt.figsize((10,5))

if (nbins_radius == 100) & (sequence == 'D305'):
    offset = 0.03
    ylim = (-0.05, 0.05)
    
if (nbins_radius == 100) & (sequence == 'A_RINGDEP_01'):
    offset = 0.02
    ylim = (-0.1, 0.1)

if (nbins_radius == 1000) & (sequence == 'D305'):
    offset = 0.1
    ylim = ((-0.1, 0.2))
    
if (nbins_radius == 1000) & (sequence == 'A_RINGDEP_01'):
    offset = 0.2
    ylim = ((-0.3, 0.3))

if (nbins_radius == 1000) & (sequence == 'D202'):
    offset = 0.2
    ylim = ((-0.3, 0.3))    
    
if (nbins_radius == 100) & (sequence == 'D202'):
    offset = 0.03
    ylim = ((-0.05, 0.05))  
    
# Plot the radial profile: mean and median

plt.plot(bins_radius / 1000, flux_mean_clean_arr, label='Mean')
plt.plot(bins_radius / 1000, flux_median_clean_arr + offset, label='Median + offset')

# Plot lines for satellite orbits

if (sequence == 'A_RINGDEP_01'):
    for name_body_i in name_body[1:]: # Loop over moons, excluding Pluto
      plt.vlines(d_pluto_body[name_body_i]/1000, -1,1, linestyle='--')
      plt.text(  d_pluto_body[name_body_i]/1000, ylim[1]*0.8, ' ' + name_body_i[0])

if (sequence == 'D305'):
    
    for rh_i in [20, 40, 60, 80]:
      plt.vlines(rh_i * r_h/1000, -1,1, linestyle='--')
      plt.text((rh_i + 2) * r_h/1000, ylim[0]*0.9, ' ' + repr(rh_i) + ' RH')
  
plt.ylim(ylim)
plt.title(sequence + ', nbins = ' + repr(nbins_radius) + \
                     ', $\delta r$ = ' + hbt.trunc(dradius,0) + ' km')

plt.xlabel('Orbital Distance [1000 km]')
plt.ylabel('DN')
plt.legend(loc = 'bottom center')

file_out = dir_out + '/profile_radial_n' + repr(nbins_radius) + '_' + sequence + '.png'
plt.savefig(file_out)
print 'Wrote: ' + file_out

plt.show()

#==============================================================================
# Make a plot of radial profile, in I/F units
#==============================================================================

ring_dn = 0.03

f_solar_1au  = 1.7e3     # erg/cm2/sec/AA, at 1 AU, at 600 nm, from http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
f_solar_pl   = f_solar / (4. * math.pi * (40.**2))

flux_median_clean_arr[0] = ring_dn # Put here the rough DN level I can measure to in cell 0. 
                                # Then we'll do all the math on that, and see what I/F this DN level converts to.

i_cgs_sr = flux_median_clean_arr / exptime / rsolar # Convert into erg/cm2/s/Angstrom/sr
i_cgs    = i_cgs_sr * pixfov  # Multiply by pixel size, to get erg/cm2/s/angstrom
i_cgs_ang= i_cgs * bandwidth # erg/cm2/sec. This is final intensity

iof = i_cgs_ang / f_solar_au

print "I/F[DN = " + repr(accuracy_dn) + "] = " + repr(iof[0])

# Concl: I/F = 1d-10. I don't believe this. It could be off by a large quantity. But it is a value.

plt.plot(bins_radius/1000, i_cgs)
plt.title(sequence)
plt.ylabel('I [cgs] [erg/cm2/s/Angstrom/sr]')
plt.xlabel('Orbital Distance [1000 km]')
plt.show()

# Now to get F, we need to divide by the solar flux
# Why is this per steradian? Becuase the ring is assumed to be resolved. 
# It's like "watts per square degree" -- a sky brightness.

# Solar flux F

# How to convert to an I/F?

stop

#==============================================================================
# Make an image showing the satellite orbits, and labled with satellite names
#==============================================================================

DO_LABEL_BODIES = False

##### A_RINGDEP_01 Sequence #####

if (sequence == 'A_RINGDEP_01'):    
    hbt.figsize((20,20))

# Make a composite image showing data, superimposed with orbits
    
    plt.subplot(1,2,1)
    
    im_composite = im.copy()
    
    mask_composite = ((mask_orbit['Charon'])   & (dist_body_pix['Charon'] > 100)) + \
                     ((mask_orbit['Styx'])     & (dist_body_pix['Styx'] > 100)) + \
                     ((mask_orbit['Hydra'])    & (dist_body_pix['Hydra'] > 100)) + \
                     ((mask_orbit['Kerberos']) & (dist_body_pix['Kerberos'] > 100)) + \
                     ((mask_orbit['Nix'])      & (dist_body_pix['Nix'] > 100))
                     
    im_composite[mask_composite] = 10
                     
    plt.imshow(stretch(im_composite))
    
    for name_body_i in name_body:
        y_pix, x_pix    = pos_body_pix[name_body_i]
        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=30)
        if (DO_LABEL_BODIES):
            plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='white', \
                     weight='bold', fontsize=15)

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    plt.title(sequence)
    
    # Plot second subplot: Un-annotated
    
    plt.subplot(1,2,2)
    plt.imshow(stretch(im_clean2))

    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)

    for name_body_i in name_body:
        y_pix, x_pix    = pos_body_pix[name_body_i]
        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=30)
    plt.title(sequence)

##### D305 Sequence #####
    
if (sequence == 'D305'):
    
    hbt.figsize((7,20))
    im_composite = im.copy()
    
    mask_composite = mask_orbit['Charon']     + mask_orbit['Hydra']        + mask_orbit['Hydra x 20'] + \
                     mask_orbit['Hydra x 40'] + mask_orbit['Hydra x 60'] + mask_orbit['Hydra x 80']
                     
    im_composite[mask_composite] = 3

    plt.subplot(1,2,1)
    plt.imshow(stretch(im_composite))
    
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    

    for name_body_i in ['Charon', 'Hydra']:
        
        y_pix, x_pix    = pos_body_pix[name_body_i]

        if (DO_LABEL_BODIES):
            plt.text(x_pix + 40, y_pix + 40, name_body_i[0], weight='bold', color='red', fontsize=12)

    plt.subplot(1,2,2)

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    plt.imshow(stretch(im_clean2))

##### D202 Sequence #####
    
if (sequence == 'D202'):
    
    hbt.figsize((7,20))
    im_composite = im.copy()
    
    mask_composite = mask_orbit['Charon']     + mask_orbit['Hydra']        + mask_orbit['Hydra x 4']


    mask_composite = ((mask_orbit['Charon'])   & (dist_body_pix['Charon'] > 50)) + \
                     ((mask_orbit['Hydra'])    & (dist_body_pix['Hydra'] > 50)) + \
                     ((mask_orbit['Hydra x 4']))

    im_composite[mask_composite] = 3

    plt.subplot(1,2,1)
    plt.imshow(stretch(im_composite))
    
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    

    for name_body_i in ['Pluto', 'Charon', 'Hydra']:
        
        y_pix, x_pix    = pos_body_pix[name_body_i]

        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=25)

        if (DO_LABEL_BODIES):
            plt.text(x_pix + 40, y_pix + 40, name_body_i[0], weight='bold', color='red', fontsize=12)

    plt.subplot(1,2,2)

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    plt.imshow(stretch(im_clean2))

plt.tight_layout()
    
file_out = dir_out + '/image_orbits_' + sequence + '.png'
plt.savefig(file_out)
print 'Wrote: ' + file_out

plt.show()

#==============================================================================
# Make a plot showing the effect of threshholding / cleaning the image
#==============================================================================

plt.set_cmap('Greys_r')

hbt.figsize((18,60))
plt.subplot(1,2,1)
plt.imshow(stretch(im))
plt.title(sequence)
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)

#plt.subplot(1,3,2)
#plt.imshow(stretch(im_clean))
#plt.title('Radius')
#plt.gca().get_xaxis().set_visible(False)
#plt.gca().get_yaxis().set_visible(False)

plt.subplot(1,2,2)
plt.imshow(stretch(im_clean2))
plt.title(sequence + ', Cleaned')
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)

file_out = dir_out + '/image_threshhold_' + sequence + '.png'
plt.savefig(file_out)
print 'Wrote: ' + file_out

plt.savefig
plt.show()

        
#==============================================================================
# XXX OLD CODE: Navigate the image using star catalog
#==============================================================================

#sequence = 'A_RINGDEP_01'
#file_raw = dir + '/ringdep_mos_v1.fits'
#
#hdulist = fits.open(file_raw)
#im = hdulist['PRIMARY'].data
#
#name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
#crval = [91, 15]  # Array with RA, Dec of center position, in degrees
#
#stars     = conesearch.conesearch(crval, 5, cache=False, catalog_db = name_cat) # center [deg], radius [deg]
#table_stars = Table(stars.array.data)
#mask = table_stars['Pmag'] < 9
#table_stars_m = table_stars[mask]  
#
#ra_stars  = table_stars_m['RAJ2000']*hbt.d2r # Convert to radians
#dec_stars = table_stars_m['DEJ2000']*hbt.d2r # Convert to radians
#
#radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
#x_stars, y_stars   = w.wcs_world2pix(radec_stars[:,0]*hbt.r2d,   radec_stars[:,1]*hbt.r2d, 0)        
#points_stars        = np.transpose((y_stars, x_stars))  # 
#
#shape = (5000,2250)
#diam_kernel=5
#
#image_1 = hbt.image_from_list_points(points_stars, shape, diam_kernel) # star catalog image: working
#image_2 = hbt.image_from_list_points(points_phot,  shape, diam_kernel) # photometry: cut off
#        
#points_phot = hbt.find_stars(im) # Weirdly, find_stars does not return magnitudes -- only positions
#            
#(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot, points_stars, np.shape(im), plot=False)


stop

# Read in one of the original, non-mosaiced images to figure out DN vs. I/F, etc.

stretch = astropy.visualization.PercentileInterval(90)

file_l2 = dir + '/mpf_0299292106_0x548_sci_2.fit'  # Level-2 from SOC
hdu_l2 = fits.open(file_l2)
im_l2  = hdu_l2['PRIMARY'].data

file_l1 = dir + '/mpf_0299292106_0x548_eng_2.fit'  # Level-1 from SOC
hdu_l1 = fits.open(file_l1)
im_l1  = hdu_l1['PRIMARY'].data

file_tod = dir + '/mpf_0299292106_02_v2.fits' # Tod gave me this. I don't know what it is.
hdu_tod = fits.open(file_tod) # )ython will not read this. "Keyword NAXIS3 not found."
im_tod = hdu_tod['PRIMARY'].data
     
im_l2 = hdu_l2['PRIMARY'].data  # 2 x 128 x 5024

# Extract some subframes to examine

im_l1_sub = im_l1[1,5:100,405:500] # Level-1. Range = 100 .. 200 mostly. Integers (ie, DN)
im_l2_sub = im_l2[1,5:100,405:500] # Level-2. Same units as l1, just scaled a bit differently. Floats (ie, calibrated DN)
rsolar    = 100190.64              # Conversion factor, DN -> cgs
 
exptime = 10

i_cgs = im_l2_sub / exptime / rsolar # Convert into erg/cm2/s/Angstrom/sr

# Comments at end of Level2 file 
#PIXSIZE =              13.0000 /Pixel size in microns                           
#READNOI =              30.0000 /Readnoise in Electrons                          
#GAIN    =              58.6000 /Gain in Electrons/DN                            
#PIXFOV  =              19.8065 /Plate scale in microrad/pix    

#PSOLAR  = '2.5541E+14'         /(DN/s)/(erg/cm^2/s/Ang), Solar spectrum. Use for unresolved
#RSOLAR  = '100190.64'          /(DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum. Use for resolved sources.
             