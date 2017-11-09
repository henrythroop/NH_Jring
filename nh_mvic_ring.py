# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 22:55:34 2016

@author: throop
"""

# This program analyzes the OUTBOUND MVIC and LORRI NH Rings data.
# (Despite this program's filename, it looks at both MVIC and LORRI.)
# This is the main program to do this analysis. It reads Tod's mosaic
# files, and generates and plots radial profiles from them.
#
# All of the values for MVIC / LORRI I/F in the data table in the Icarus Pluto Ring paper
# are directly from this code. In order to get these values, set the file (e.g., MVIC 202) and binning (e.g., 1000)
# and run.
#
# HBT Nov 2016 thru May 2017.

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
import matplotlib # So I can use matplotlib.rc
import numpy as np
import astropy.modeling
import astropy.visualization
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import astropy.units as u
import astropy.constants as c

import spiceypy as sp
import skimage
#from   itertools   import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy     import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils   import daofind
#import wcsaxes
import time
from scipy.interpolate import griddata

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

import hbt

dir_mvic      = '/Users/throop/Data/NH_MVIC_Ring'    # MVIC
dir_lorri     = '/Users/throop/Data/NH_LORRI_Ring'   # LORRI


file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

sp.furnsh(file_tm) # Start up SPICE

#==============================================================================
# Initialize constants
#==============================================================================

# '-new-image' is added by astrometry.net
# '_fixed'     is added by me to indicate that I have fixed the header (that it, added SCET field)
#files = ['mvic_d305_sum_mos_v1.fits', 'ringdep_mos_v1.fits']

#file = dir + '/mvic_d305_sum_mos_v1.fits'

# SELECT THE MVIC OR LORRI MOSAIC TO ANALYZE

#sequence        = 'D305' # MVIC
#sequence        = 'O_RINGDEP_A_1'  # Departure imaging, closest MVIC image of the whole system
#sequence        = 'D202' # MVIC
#sequence        = 'D211'  # MVIC
sequence        = 'U_TBD_4'   # MVIC

#sequence        = 'D202_LORRI'
#sequence        = 'D305_LORRI'

# SELECT THE NUMBER OF RADIAL BINS TO USE

nbins_radius = 100
#nbins_radius = 1000
 
DO_FIX_FITS     = False
DO_ANALYZE      = True

DO_PLOT_TITLE  = False   # Plot a title on the radial profile. Turn off for publication.
#DO_PLOT_TITLE = True

#PIXSIZE =              13.0000 /Pixel size in microns                           
#READNOI =              30.0000 /Readnoise in Electrons                          
#GAIN    =              58.6000 /Gain in Electrons/DN                            
#PIXFOV  =              19.8065 /Plate scale in microrad/pix    

#PSOLAR  = '2.5541E+14'         /(DN/s)/(erg/cm^2/s/Ang), Solar spectrum. Use for unresolved
#RSOLAR  = '100190.64'          /(DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum. Use for resolved sources.

if ('LORRI' in sequence):
    IS_LORRI, IS_MVIC = True, False

else:
    IS_MVIC, IS_LORRI = (True, False)

# Set the directory based on which instrument we are using

if (IS_LORRI):
    dir = dir_lorri
else:
    dir = dir_mvic    

dir_out = dir + '/out'

fontsize = 15    # Font size for radial profile plots
figsize = (10,5) # Figure size for radial profile plots

r_pluto_km = 1187

m_pluto = 1.309e22 * u.kg
r_hill = 33*u.AU * (m_pluto / (3 * c.M_sun))**(1/3)

#==============================================================================
# Initialize parameters for each of the possible mosaics we can analyze
#==============================================================================

# *_wcs.fits file are created by astrometry.net. Their default name is new-image.fits .
# *_header.fits indicates that I have added ET and any other necessary header info into the file. 
# *_pl.fits indicates that I have added backplanes.

if (sequence == 'D305'):
  file_wcs    = dir + '/mvic_d305_sum_mos_v1_wcs.fits' # Load the navigated image, with WCS  
  utc = '2015::305 00:00:00'  # Set the time 
  stretch = astropy.visualization.PercentileInterval(99)
  
if (sequence == 'O_RINGDEP_A_1'):
  file_wcs = dir + '/ringdep_mos_v1_wcs.fits'
  utc = '2015 Jul 15 18:50:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)

if (sequence == 'D202'):
  file_wcs = dir + '/mvic_d202_mos_v1_wcs.fits'
  utc = '2015::202 00:00:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
  
if (sequence == 'D211'):
  file_wcs = dir + '/mvic_d211_mos_v1_wcs.fits'
  utc = '2015::211 00:00:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
  
if (sequence == 'D202_LORRI'):
  file_wcs = dir + '/ring_lor202_mos_v1_wcs.fits'
  utc = '2015::202 00:00:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
    
if (sequence == 'D305_LORRI'):
  file_wcs = dir + '/ring_lor305_mos_v3_wcs.fits'
  utc = '2015::305 00:00:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
  

file_header    = file_wcs.replace('.fits', '_header.fits')    # File with fixed FITS ET info
file_header_pl = file_header.replace('.fits', '_pl.fits')     # File for backplanes

#==============================================================================
# Repair / edit the FITS files if needed
#==============================================================================
    
if (DO_FIX_FITS):

    file_image_in = file_wcs
   
    file_header_out = file_header
    
# Add a missing header field and write out

    if (sequence == 'O_RINGDEP_A_1'):
        file_header_in = dir + '/mpf_0299292106_0x548_sci_2.fit'
    if (sequence == 'D202'):
        file_header_in = dir + '/mpf_0299793106_0x539_sci_2.fit'
    if (sequence == 'D211'):
        file_header_in = dir + '/mpf_0300584506_0x539_sci_3.fit'
    if (sequence == 'D305'):
        file_header_in = dir + '/mpf_0308706706_0x539_sci_6.fit'

    print("Reading file: " + file_header_in)
    print("Reading file: " + file_image_in)
    
    hdulist_header_in = fits.open(file_header_in)
    hdulist_image_in  = fits.open(file_image_in) # Has Tod's image, and WCS data

# Start with Tod's mosaics, which have the correct image and the correct WCS.
# Append all of the regular NH FITS info to the end of that, from a raw SOC image.
    
    for card in hdulist_header_in['PRIMARY'].header.cards:
        hdulist_image_in['PRIMARY'].header.append(card)

# Now write out the image data (ie, Tod's mosaics)
    try:
        hdulist_image_in['PRIMARY'].header.remove('NAXIS3') # Remove a spurious keyword that Tod put there
    except ValueError:
        pass
    
    hdulist_image_in.writeto(file_header_out, clobber=True)
    
    print("Wrote new FITS file with fixed header: " + file_header_out)

    hdulist_header_in.close()
    hdulist_image_in.close()
    
# Create the backplanes. 
# Start by reading in the file we just created.

    print("Generating backplanes...")
    
    planes = hbt.create_backplane(file_header_out, frame = 'IAU_PLUTO', 
                                  name_target='Pluto', name_observer='New Horizons')

# Create backplaned FITS file.

# Create a new FITS file, made of an existing file plus these new planes

    hdulist = fits.open(file_header_out)  # Read in the file we've just created, which has full fixed header

    type_out = 'float32'  # Write backplane in single precision, to save space.
    hdu_out = fits.HDUList() # Create a brand new FITS file
    hdu_out.append(hdulist['PRIMARY'])
    hdu_out.append(fits.ImageHDU(planes['RA'].astype(type_out), name = 'RA'))
    hdu_out.append(fits.ImageHDU(planes['Dec'].astype(type_out), name = 'Dec'))
    hdu_out.append(fits.ImageHDU(planes['Phase'].astype(type_out), name = 'Phase'))
    hdu_out.append(fits.ImageHDU(planes['Longitude_eq'].astype(type_out), name = 'Longitude_eq'))
    hdu_out.append(fits.ImageHDU(planes['Radius_eq'].astype(type_out), name = 'Radius_eq'))
    
    hdu_out.writeto(file_header_pl, clobber=True)
    print("Wrote file: " + file_header_pl)

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
utc = sp.et2utc(et, 'C', 0)

w = WCS(file)

# Calculate the sub-observer latitude (ie, ring tilt angle)

(vec, lt) = sp.spkezr('New Horizons', et, 'IAU_PLUTO', 'LT', 'Pluto')
vec = vec[0:3]
(junk, lon_obs, lat_obs) = sp.reclat(vec) # Get latitude, in radians

#==============================================================================
# Make a plot of the image + backplanes
#==============================================================================

plt.set_cmap('Greys_r')

hbt.figsize((9,6)) # This gets ignored sometimes here. I'm not sure why??
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
# Make a labled plot with Pluto, Charon, Nix, Hydra, etc.
#==============================================================================

hbt.figsize((5,12))

offset_x, offset_y = 100,100  # Offset PCNHSK labels by this many pixels

color = 'red'

if (sequence == 'D305_LORRI'):
    
    offset_x, offset_y = 20,20
    hbt.figsize((12,15))
    color='white'

if (sequence == 'D202_LORRI'):
    color = 'white'
    hbt.figsize((12,15))
    offset_x, offset_y = 20,20    
    
plt.imshow(stretch(im))
    
pos_body_pix = {}  # Create a new dictionary to save positions of each body

name_body = ['Pluto', 'Charon', 'Nix', 'Hydra', 'Styx', 'Kerberos']
for name_body_i in name_body:
    vec,lt = sp.spkezr(name_body_i, et, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    (junk,ra,dec) = sp.recrad(vec_sc_targ) # Get the RA / Dec of the object
    x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
    plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
             markeredgecolor=color, markeredgewidth=2)
    plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color=color, fontsize=15)

    pos_body_pix[name_body_i] = np.array([y_pix, x_pix])  # Save position of each body as a dictionary entry tuple
        
plt.title(sequence)
plt.xlim((0,np.shape(im)[1]))
plt.ylim((0,np.shape(im)[0]))

plt.tight_layout()
file_out = dir_out + '/image_annotated_' + sequence + '.png'
plt.savefig(file_out)
print('Wrote: ' + file_out)

plt.show()

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

dn_mean_arr   = np.zeros(nbins_radius)
dn_median_arr = np.zeros(nbins_radius)
dn_std_arr    = np.zeros(nbins_radius)
npix_arr        = np.zeros(nbins_radius)

dn_mean_clean_arr   = np.zeros(nbins_radius)
dn_median_clean_arr = np.zeros(nbins_radius)

for i in range(nbins_radius-1):
    is_good = np.array(radius > bins_radius[i]) \
            & np.array(radius < bins_radius[i+1]) \
            & (dist_body_pix['Pluto']  > 100) \
            & (dist_body_pix['Charon'] > 100) \
            & (im != 0) \
            & (im != -10.)
            
            # Mask out the zero pixels (e.g., between MVIC mosaic frames) 
            # Also need to mask out area near Pluto and near Charon
            
    bin_number_2d[is_good] = i          # For each pixel in the array, assign a bin number

    dn_mean_arr[i]         = np.mean(im[is_good])
    dn_mean_clean_arr[i]   = np.mean(im_clean2[is_good])
    dn_median_arr[i]       = np.median(im[is_good])
    dn_median_clean_arr[i] = np.median(im_clean2[is_good])
    dn_std_arr[i]          = np.std(im[is_good]) # stdev
    
    npix_arr[i] = np.sum(is_good) # Count of number of pixels in this radial bin

# Now normalize these so they have median at zero. Could argue whether to do this or not, but I think it's fine.

dn_mean_arr         -= np.nanmedian(dn_mean_arr)
dn_mean_clean_arr   -= np.nanmedian(dn_mean_clean_arr)
dn_median_arr       -= np.nanmedian(dn_mean_clean_arr)
dn_median_clean_arr -= np.nanmedian(dn_median_clean_arr)

    
#==============================================================================
# Calc orbital distances (from Pluto) for each body
#==============================================================================

d_pluto_body = {}
mask_orbit = {}

if (sequence == 'D305'):
    d_radius = 4000  # Halfwidth to plot, of orbit, in km
  
if (sequence == 'O_RINGDEP_A_1'):
    d_radius = 100
  
if (sequence == 'D202'):
    d_radius = 500

if (sequence == 'D211'):
    d_radius = 1000

if (sequence == 'D202_LORRI'):
    d_radius = 500

if (sequence == 'D305_LORRI'):
    d_radius = 4000
    
# Look up the distance using SPICE

for name_body_i in name_body:
  (vec_body,junk) = sp.spkezr(name_body_i, et, 'IAU_PLUTO', 'LT', 'Pluto')
  d_pluto_body[name_body_i] = sp.vnorm(vec_body[0:3])
   
# Generate a pixel mask showing the orbit of each body

  mask_orbit[name_body_i] = \
    np.array(radius > (d_pluto_body[name_body_i] - d_radius)) & \
    np.array(radius < (d_pluto_body[name_body_i] + d_radius))

    
r_h = d_pluto_body['Hydra']  # Hydra orbital radius

mask_orbit['Hydra x 2'] = \
    np.array(radius > (r_h*2 - d_radius)) & np.array(radius < (r_h*2 + d_radius))
    
mask_orbit['Hydra x 4'] = \
    np.array(radius > (r_h*4 - d_radius)) & np.array(radius < (r_h*4 + d_radius))
    
mask_orbit['Hydra x 5'] = \
    np.array(radius > (r_h*5 - d_radius)) & np.array(radius < (r_h*5 + d_radius))
    
mask_orbit['Hydra x 6'] = \
    np.array(radius > (r_h*6 - d_radius)) & np.array(radius < (r_h*6 + d_radius))
        
mask_orbit['Hydra x 10'] = \
    np.array(radius > (r_h*10 - d_radius)) & np.array(radius < (r_h*10 + d_radius))

mask_orbit['Hydra x 20'] = \
    np.array(radius > (r_h*20 - d_radius)) & np.array(radius < (r_h*20 + d_radius))

mask_orbit['Hydra x 30'] = \
    np.array(radius > (r_h*30 - d_radius)) & np.array(radius < (r_h*30 + d_radius))

mask_orbit['Hydra x 40'] = \
    np.array(radius > (r_h*40 - d_radius)) & np.array(radius < (r_h*40 + d_radius))

mask_orbit['Hydra x 60'] = \
    np.array(radius > (r_h*60 - d_radius)) & np.array(radius < (r_h*60 + d_radius))
    
mask_orbit['Hydra x 80'] = \
    np.array(radius > (r_h*80 - d_radius)) & np.array(radius < (r_h*80 + d_radius))    

#==============================================================================
# Convert measurements from DN into I/F and optical depth
#==============================================================================

dist_au    = hdulist[0].header['SPCTSORN']*u.km.to('AU')  # Heliocentric distance in AU. Roughly 32 AU.

lat_subsc  = hdulist[0].header['SPCTSCLA']                 # Sub-sc Lat, in deg. Roughly 43 deg.

RSOLAR_l11 = 221999.98             # (DN/s/pixel)/(erg/cm^2/s/A/sr) for LORRI 1x1
#RSOLAR_l44_OLD = 3800640.0             # (DN/s/pixel)/(erg/cm^2/s/A/sr) for LORRI 4x4. = 221999.98 * 1.07 * 16.
RSOLAR_l44 = 4.09e6                # Revised LORRI 4x4 value, 22-May-2017, Hal Weaver. 
#RSOLAR_mpf_OLD = 100190.64             # (DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum. MVIC Pan Frame. 
RSOLAR_mpf = 9.8813e4
FSOLAR_LORRI  = 176.	     	    # Solar flux at r=1 AU in ergs/cm^2/s/A
FSOLAR_MPF    = 145.                # Solar flux, assuming pivot wavelength 6920 A -- see Hal email 18-May-2017  

n_sig      = 3                     # Do we do 3sigma? or 4 sigma? Set it here (3 or 4)

if IS_LORRI:
    exptime = 0.4  # Tod's LORRI mosaics have 200ms and 400ms in them, and are normalized to 400 ms.

    FSOLAR  = FSOLAR_LORRI
    
    if hdulist[0].header['SFORMAT'] == '1X1':
        RSOLAR = RSOLAR_l11
        
    if hdulist[0].header['SFORMAT'] == '4X4':
        RSOLAR = RSOLAR_l44
                
if IS_MVIC:
    exptime = 10   # Tod's MVIC Pan Frame mosaics are all 10 sec exposures, and averaged (not summed).
    RSOLAR = RSOLAR_mpf
    FSOLAR = FSOLAR_MPF

#==============================================================================
# Dummy check: Calculate sub-sc latitude using SPICE, rather than from FITS header [concl: working fine]
#==============================================================================

(vec,ltime) = sp.spkezr('Pluto', et, 'J2000', 'LT+S', 'New Horizons')
sp.illum('Pluto', et, 'LT+S', 'New Horizons')
(subpnt_rec, junk, junk) = sp.subpnt('Intercept: ellipsoid', 'Pluto', et, 'IAU_PLUTO', 'LT+S', 'New Horizons')
(rad_spice, lon_subsc_spice, lat_subsc_spice) = sp.reclat(subpnt_rec)

print("Sub-sc lat = {} [FITS]; {} [SPICE]".format(lat_subsc, hbt.r2d * lat_subsc_spice))

#==============================================================================
# Convert from DN to irradiance (ergs/cm2/s/sr/A). From Hal's writeup.
#==============================================================================
 
dn = dn_median_clean_arr
   
I = dn / exptime / RSOLAR   

# Convert from irradiance to I/F. From Hal's writeup.

iof = math.pi * I * dist_au**2 / FSOLAR

# Convert from I/F to 'normal I/F'.
# This value -- 4 mu I/F -- is exactly the same as tau \varpi P, according to Throop 2004 eq 1

# Get ring emission angle. This should use 90-lat_subsc, rather than lat_subsc that I used in original submitted paper.
# I should use the emission angle... that is, from *normal*, not from *equator*.
# Error discovered 29-Aug-2017. Fixed 21-Sep-2017.
# Correct usage is given by TPW04 @ 63: 
#   I/F = tau pi_0 P / (4 mu)
#   mu  = |cos(e)|
#   e   = emission angle. 0 for directly above ring.
#   B   = elevation angle = 90 - e.  [Same as lat_subsc here.]
#         So, mu = cos(e) = cos(90-B) = cos(90-lat_subsc)    

mu = math.cos(math.pi/2 - lat_subsc * hbt.d2r)  # mu = cos(lat)
iof_normal = 4 * mu * iof  # (I/F)_normal = 4 mu I/F

# Calculate the 3sigma iof value. For this, should we use the inner data points only (yes, most of the time)
# or take stdev of all data points (which will include some extrema which I think should be tossed.)
# Only do this for the MVIC frames. The LORRI frames have bad stray light throughout, and we want to take the stdev
# across the whole frame, not just the inner region.

DO_STDEV_INNER_ONLY = True

#DO_STDEV_INNER_ONLY = False

if (DO_STDEV_INNER_ONLY) and (IS_MVIC):
     x0 = int(np.size(iof_normal)*0.25)
     x1 = int(np.size(iof_normal)*0.75)
    
else:
    x0 = 0
    x1 = hbt.sizex(iof_normal)
     
std_iof_normal = np.nanstd(iof_normal[x0:x1])   

# Special case: for the O_RINGDEP_A_1 sequence, the edges are bad, but the signal 
# in the middle is good. So, when we take the SNR, do it *only* of the middle region.
# We will show the full plot and it will be obvious what we did.

#if (sequence == 'O_RINGDEP_A_1'):
#     x0 = int(np.size(iof_normal)*0.25)
#     x1 = int(np.size(iof_normal)*0.75)
#     std_iof_normal = np.nanstd(iof_normal[x0:x1])   
    
#==============================================================================
# Make plots of radial profile
#==============================================================================

matplotlib.rc('font', size = fontsize)

hbt.figsize((10,5))
    
if (nbins_radius == 100) & (sequence == 'O_RINGDEP_A_1'):
    offset = 0.02
    ylim_dn = (-0.2, 0.2)

if (nbins_radius == 1000) & (sequence == 'O_RINGDEP_A_1'):
    offset = 0.
    ylim_dn = ((-0.2, 0.2))

if (nbins_radius == 100) & (sequence == 'D202'):
    offset = 0.2
    ylim_dn = ((-0.15, 0.15))    
    
if (nbins_radius == 1000) & (sequence == 'D202'):
    offset = 0.2
    ylim_dn = ((-0.15, 0.15))    
    
if (nbins_radius == 100) & (sequence == 'D211'):
    offset = 0.03
    ylim_dn = ((-0.05, 0.05))  
        
if (nbins_radius == 1000) & (sequence == 'D211'):
    offset = 0.06
    ylim_dn = ((-0.05, 0.05))  
        
if (nbins_radius == 100) & (sequence == 'D305'):
    offset = 0.03
    ylim_dn = (-0.05, 0.05)

if (nbins_radius == 1000) & (sequence == 'D305'):
    offset = 0.1
    ylim_dn = ((-0.1, 0.2))

if (nbins_radius == 1000) & (sequence == 'D202_LORRI'):
    offset = 0.1
    ylim_dn = ((-1, 1))

if (nbins_radius == 100) & (sequence == 'D202_LORRI'):
    offset = 0.1
    ylim_dn = ((-1, 1))
    
if (nbins_radius == 1000) & (sequence == 'D305_LORRI'):
    offset = 0.1
    ylim_dn = ((-2, 2))

if (nbins_radius == 100) & (sequence == 'D305_LORRI'):
    offset = 0.1
    ylim_dn = ((-2, 2))

# Plot the radial profile: mean and median

DO_PLOT_PROFILE_MEDIAN = True
DO_PLOT_PROFILE_MEAN   = False

hbt.figsize((10,6))

fig, ax1 = plt.subplots()

if (DO_PLOT_PROFILE_MEAN):
    ax1.plot(bins_radius / 1000, dn_mean_clean_arr, label='Mean')
    
if (DO_PLOT_PROFILE_MEDIAN):

    ax1.plot(bins_radius / 1000, 
             dn_median_clean_arr +
             offset*DO_PLOT_PROFILE_MEAN, # Ignore offset if plotting only one curve
             label='Median')

# Plot lines for satellite orbits

alpha = 0.3

if (sequence == 'O_RINGDEP_A_1'):
    for rh_i in [1]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 0.05) * r_h/1000, ylim_dn[0]*0.9, ' ' + repr(rh_i) + ' R$_H$')
      
if (sequence == 'D305'):
    for rh_i in [40, 80]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 2) * r_h/1000, ylim_dn[0]*0.9, ' ' + repr(rh_i) + ' R$_H$')

if (sequence == 'D211'):
    for rh_i in [5, 10]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 0.5) * r_h/1000, ylim_dn[0]*0.9, ' ' + repr(rh_i) + ' R$_H$')
      
if (sequence == 'D202'):
    for rh_i in [1, 2, 4]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 0.1) * r_h/1000, ylim_dn[0]*0.8, ' ' + repr(rh_i) + ' $R_H$')

if (sequence == 'D202_LORRI'):
    for rh_i in [1, 2]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 0.05) * r_h/1000, ylim_dn[0]*0.8, ' ' + repr(rh_i) + ' $R_H$')
      
if (sequence == 'D305_LORRI'):
    for rh_i in [5, 10, 20, 40]:
      plt.vlines(rh_i * r_h/1000, -10,10, linestyle='--', alpha=alpha)
      plt.text((rh_i + 0.1) * r_h/1000, ylim_dn[0]*0.8, ' ' + repr(rh_i) + ' $R_H$')
            
ax1.set_ylim(ylim_dn)

# Create a second y axis

# We don't plot a second time. We just set the y range on the second axis.
# This line below looks funny, but it calculates the exact ratio of dn:iof, and sets ylimit based on that.

ylim_iof = np.array(ylim_dn) * np.nanmedian(iof_normal / dn)

ax2 = ax1.twinx() # Yes, twinx means they will share x axis (and have indep y axes) - correct as written.
#ax2.plot(bins_radius/1000, iof_normal)

# Set scaling factor for I/F

if (IS_LORRI):
    iof_units = 1e-5
    iof_units_str = r"$\times 10^{-5}$"
    
if (IS_MVIC):
    iof_units = 1e-6    
    iof_units_str = r"$\times 10^{-6}$"

ax2.set_ylim(ylim_iof / iof_units)
ax2.set_ylabel('Normal I/F [{}]'.format(iof_units_str))

if (DO_PLOT_TITLE):
    ax1.set_title(sequence + ', nbins = ' + repr(nbins_radius) + \
#                     ', $\delta r$ = {:.0f} km'.format(dradius,0) +                     
                     ', {}$\sigma$ I/F $\leq$ {:.1e}'.format(n_sig, n_sig*std_iof_normal))

ax1.set_xlabel('Orbital Distance [1000 km]')
ax1.set_ylabel('DN')

# Plot the legend, iff we are plotting two lines (mean and median)

if (DO_PLOT_PROFILE_MEDIAN and DO_PLOT_PROFILE_MEAN):
    ax1.legend(loc = 'lower center') # also 'lower center'

# Plot a text with the final derived I/F limit
#ax1.text(0.02, 0.8, '3$\sigma$ I/F = {:.1e}'.format(3*std_iof_normal), transform=ax1.transAxes)

file_out = dir_out + '/profile_radial_n' + repr(nbins_radius) + '_' + sequence + \
    ('' if DO_PLOT_TITLE else '_notitle') + '.png'

fig.savefig(file_out)
print('Wrote: ' + file_out)

fig.show()

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
    
##### D305 LORRI Sequence #####
    
if (sequence == 'D305_LORRI'):
    
    hbt.figsize((7,20))
    im_composite = im.copy()
    
    mask_composite = mask_orbit['Charon']     + mask_orbit['Hydra']        + mask_orbit['Hydra x 5'] + \
                     mask_orbit['Hydra x 10'] + mask_orbit['Hydra x 20'] + mask_orbit['Hydra x 30']
                     
    im_composite[mask_composite] = 3

#    im_composite[im_composite == im_composite[0,0]] = np.amin(im_composite)
   
    plt.subplots() 
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
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))

    plt.tight_layout()
    
##### D202 LORRI Sequence #####
    
if (sequence == 'D202_LORRI'):
    
    hbt.figsize((20,20))
    im_composite = im.copy()
    
#    mask_composite = mask_orbit['Charon']     + mask_orbit['Hydra']        + mask_orbit['Hydra'] + \
#                     mask_orbit['Hydra x 2']

    mask_composite = ((mask_orbit['Charon'])   & (dist_body_pix['Charon'] > 20)) + \
                     ((mask_orbit['Styx'])     & (dist_body_pix['Styx'] > 20)) + \
                     ((mask_orbit['Hydra'])    & (dist_body_pix['Hydra'] > 20)) + \
                     ((mask_orbit['Kerberos']) & (dist_body_pix['Kerberos'] > 20)) + \
                     ((mask_orbit['Nix'])      & (dist_body_pix['Nix'] > 20)) + \
                     mask_orbit['Hydra x 2']
                     
                     
#    im_composite[mask_composite] = 10
                     
    im_composite[mask_composite] = 10

#    im_composite[im_composite == im_composite[0,0]] = np.amin(im_composite)
    
    plt.subplots()
    plt.subplot(1,2,1)
    plt.imshow(stretch(im_composite))
    
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)

    for name_body_i in name_body:
        y_pix, x_pix    = pos_body_pix[name_body_i]
        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=30)
        if (DO_LABEL_BODIES):
            plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='white', \
                     weight='bold', fontsize=15)

            
#    for name_body_i in ['Charon', 'Hydra']:
#        
#        y_pix, x_pix    = pos_body_pix[name_body_i]
#
#        if (DO_LABEL_BODIES):
#            plt.text(x_pix + 40, y_pix + 40, name_body_i[0], weight='bold', color='red', fontsize=12)
#
    plt.subplot(1,2,2)

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    plt.imshow(stretch(im_clean2))
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    
    plt.tight_layout()
        
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

##### D211 Sequence #####
    
if (sequence == 'D211'):
    
    DO_LABEL_BODIES = False
    
    hbt.figsize((7,20))
    im_composite = im.copy()

    mask_composite = ((mask_orbit['Hydra'])    & (dist_body_pix['Hydra'] > 50)) + \
                     ((mask_orbit['Hydra x 10'])) + \
                     ((mask_orbit['Hydra x 5']))

    im_composite[mask_composite] = 3

    plt.subplot(1,2,1)
    plt.imshow(stretch(im_composite))
    
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    

    for name_body_i in ['Pluto', 'Hydra']:
        
        y_pix, x_pix    = pos_body_pix[name_body_i]

        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=25)

        if (DO_LABEL_BODIES):
            plt.text(x_pix + 40, y_pix + 40, name_body_i[0], weight='bold', color='red', fontsize=12)

    plt.subplot(1,2,2)

    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.title(sequence)
    
    plt.imshow(stretch(im_clean2))
    
#plt.tight_layout()
    
file_out = dir_out + '/image_orbits_' + sequence + '.png'
plt.savefig(file_out)
print('Wrote: ' + file_out)

plt.show()

#==============================================================================
# Make a plot showing the effect of threshholding / cleaning the image
#==============================================================================

plt.set_cmap('Greys_r')

if (IS_MVIC):
    hbt.figsize((18,60))

if (IS_LORRI):
    hbt.figsize((25,35))
    
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

plt.tight_layout()

file_out = dir_out + '/image_threshhold_' + sequence + '.png'
plt.savefig(file_out)
print('Wrote: ' + file_out)

plt.savefig
plt.show()

#==============================================================================
# Make a data table about this image and sequence
#==============================================================================

r_hill_km = r_hill.to('km').value

#nsig = 3   # 3 sigma, 1 sigma, etc.

print("Sequence = {}".format(sequence))
print("N = {} bins".format(nbins_radius))
print("Radius = {:.1f} .. {:.1f} km".format(np.amin(radius), np.amax(radius)))
print("Radius = {:.1f} .. {:.1f} r_pluto".format(np.amin(radius)/r_pluto_km, np.amax(radius)/r_pluto_km))
print("Radius = {:.1f} .. {:.1f} R_hydra".format(np.amin(radius)/r_h, np.amax(radius)/r_h))
print("Radius = {:.1f} .. {:.2f} R_hill".format(np.amin(radius)/r_hill_km, np.amax(radius)/r_hill_km))
print("Radial resolution = {:.1f} km".format( (np.amax(radius)-np.amin(radius)) / nbins_radius  ))
print("{} sigma I/F = {:.2e}".format(n_sig, n_sig*std_iof_normal))
