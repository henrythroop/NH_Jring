#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 22:45:55 2016

@author: throop
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy
from   astropy.table import Table
import astropy.table as table
from   astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import math
import hbt
import spiceypy as sp
import astropy.units as u
import astropy.constants as c
import pickle
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
import astroquery
from   astroquery.vizier import Vizier 
from matplotlib.patches import Ellipse
import os.path # For expanduser


# NB: In future, conesearch will move from astropy to astroquery -- see mailing list 17-Mar-2017.
# But it hasn't happened yet, so even if I wanted to be pro-active, I can't.

#==============================================================================
# A quick wrapper to query the Gaia catalog using Vizier
# From https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
#==============================================================================

def gaia_query(ra_deg, dec_deg, rad_deg, maxmag=20, 
               maxsources=10000): 
    """
    Query Gaia DR1 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field 
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'phot_g_mean_mag'], 
                    column_filters={"phot_g_mean_mag": 
                                    ("<%f" % maxmag)}, 
                    row_limit = maxsources) 
 
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, 
                           unit=(u.deg, u.deg), 
                           frame='icrs')
    return vquery.query_region(field, 
                               width=("%fd" % rad_deg), 
                               catalog="I/337/gaia")[0] 
    

#==============================================================================
# Return a list of NH MegaCam Gaia stars that match a given position and mag limit
#==============================================================================
# NB: After I finished this function, I realized that RA of catalog excerpt is
# different than RA of MU69 KEM encounter. So maybe this section is just academic.

def nh_gaia_query(ra_deg, dec_ddeg, rad_deg, maxmag=20, maxsources=10000):
    
    """
    Return a list of NH MegaCam Gaia stars that match a given position and mag limit
    """
    
    # Load the catalog into an array.
    # In an ideal world we would not load this if it was already loaded.
    # But I don't know how to do that. Function only gets called 1x per execution, so 
    # maybe not a big deal?
    
    nh_gaia = load_nh_gaia()
    
    # Do a set of boolean matches to extract and return the proper stars.
    # This logic here is simplistic: square box, ignore cos(dec), ignore 360->0 crossing, etc.
    # But for the case of NH KEM, that is sufficient.
    
    is_ra = np.logical_and( nh_gaia['RA'] < (ra_deg + rad_deg),
                            nh_gaia['RA'] > (ra_deg - rad_deg) )

    is_dec = np.logical_and(nh_gaia['Dec'] < (dec_deg + rad_deg),
                            nh_gaia['Dec'] > (dec_deg - rad_deg) )
    
    is_mag_g = np.logical_and( nh_gaia['g'] < maxmag,
                               nh_gaia['g'] > 0)
    
    is_mag_r = np.logical_and( nh_gaia['r'] < maxmag,
                               nh_gaia['r'] > 0)

    # Combine all these bitmasks
    
    is_good = np.logical_and(
                  np.logical_or(is_mag_g, is_mag_r),
                  np.logical_and(is_ra, is_dec) )


    # Return the selected values
    
    return nh_gaia[is_good]

#==============================================================================
# Read the NH MegaCam-Gaia star catalog from disk. Use a pickle file if available.
#==============================================================================
# NB: After I finished this function, I realized that RA of catalog excerpt is
# different than RA of MU69 KEM encounter. So maybe this section is just academic.
    
def load_nh_gaia():

    """
    Read the NH MegaCam-Gaia star catalog from disk. 
    Use a pickle file if available. If not, read text file 
    and then create a pickle file.
    """
    
    dir_gaia = os.path.expanduser('~') + '/Data/Catalogs/Gaia/' 
    file_txt = 'nh16.mega.gaia.rdmpm.txt'
    file_pickle = file_txt.replace('.txt', '.pkl')
    
    if os.path.isfile(dir_gaia + file_pickle):
        lun = open(dir_gaia + file_pickle, 'rb')
        gaia = pickle.load(lun)
        lun.close()
        print("Loaded: " + file_pickle)
    
    else:
        
        # Read the NH MegaCam-Gaia catalog from disk
      
        gaia = astropy.io.ascii.read(dir_gaia + file_txt, format = 'basic')
          
        # Make a plot of the Gaia stars
          
        plt.plot(gaia['RA'], gaia['Dec'], linestyle='none', marker = '.', ms=0.005)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title('NH MegaCam-Gaia catalog')
        plt.show()
      
      # Save it as a pickle file
      
        lun = open(dir_gaia + file_pickle, 'wb')
        pickle.dump(gaia, lun) 
        print("Wrote: " + dir_gaia + file_pickle)
        lun.close()
    
    return gaia
  
#==============================================================================
# Define kernels and files
#==============================================================================

file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69.tm'       # C/A time 11:00:00 (??) - older 
file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69_tcm22.tm' # C/A time 07:00:00 - newer version

file_hd_pickle = '/Users/throop/git/NH_rings/cat_hd.pkl'

#==============================================================================
# Initialize settings
#==============================================================================

sp.furnsh(file_tm) # Start up SPICE

hour      = 3600
day       = 24 * hour

fs        = 15   # General font size

if ('tcm22' in file_tm):
    utc_ca    = '2019 1 jan 07:00:00'
else:
    utc_ca    = '2019 1 jan 03:09:00' # MU69 C/A time. I got this from GV.
    
et_ca     = sp.utc2et(utc_ca)
dt_tof    = 180                   # Time-of-flight uncertainty (halfwidth, seconds)

hbt.figsize((15,10))
mag_limit = 8 # Plot stars brighter than this.
plot_tick_every = 360  # Plot a time-tick every __ seconds

DO_PLOT_TOF = False               # In general, plotting TOF is not necessary -- it
                                  # is a translation in position along track,
                                  # and changes the occultation times but not positions.

DO_PLOT_LORRI = True              # Plot limits of LORRI FOV?
          
DO_SCALEBAR_ARCSEC = False

DO_UNITS_TITLE_DAYS = False       # Flag: Use Days or Hours for the units in the plot title  
        
DO_UNITS_TICKS_DAYS = False

DO_PLOT_UNCERTAINTY_MU69 = False  # Make a plot of errorbars in MU69 position?
          
# Define the start and end time for the plot. This is the main control to use.

case = 6

if (case == 1): # Outbound, 0.3 .. 96 hour
    et_start  = et_ca + 0.3*hour
    et_end    = et_ca + 96*hour
    pad_ra_deg  = 1 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 1
    DO_LABEL_HD = True # Plot star IDs on chart
    DO_PLOT_HD = True
    hbt.figsize((15,10))
    mag_limit = 8 # Plot stars brighter than this.

if (case == 2): # Inbound, -96h .. -0.1h
    et_start  = et_ca -96*hour
    et_end    = et_ca -0.1*hour
    pad_ra_deg  = 2 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 2
    DO_PLOT_HD  = True
    DO_LABEL_HD = True
    hbt.figsize((25,10))
    mag_limit = 8 # Plot stars brighter than this.

if (case == 3): # Outbound, 2 .. 200 hour
    et_start  = et_ca + 2*hour
    et_end    = et_ca + 200*hour
    pad_ra_deg  = 0.5 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.5
    DO_PLOT_HD = True
    DO_LABEL_HD = True # Plot star IDs on chart
    hbt.figsize((15,10))
    mag_limit = 10 # Plot stars brighter than this.

if (case == 4): # Outbound, 2d .. 10 *** This one is key!! It has one good occultation.
    et_start  = et_ca + 2*day
    et_end    = et_ca + 20*day
    pad_ra_deg  = 0.01 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.01
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False

if (case == 4.5): # Outbound, zoom on the good one!
    et_start  = et_ca + 3.0*day
    et_end    = et_ca + 20*day
    pad_ra_deg  = 0.005 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.005
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False
    
if (case == 5): # Inbound, 2d .. 10
    et_start  = et_ca - 20*day
    et_end    = et_ca - 2*day
    pad_ra_deg  = 0.01 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.01
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False

    
if (case == 6): # Inbound, K-40d .. K-14 (for Spencer one-off project, not for an MT)
                # Using K-40 since it is one that JS did calculations for in his 22-Mar-2017 email.
                
    et_start         = et_ca - 40*day
    et_end           = et_ca - 14*day
    pad_ra_deg       = 0.008 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg      = 0.008
    
    DO_PLOT_HD       = False
    DO_LABEL_HD      = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO     = True
    DO_LABEL_USNO    = True # DO_PLOT_USNO
    DO_PLOT_GAIA     = True
    DO_LABEL_GAIA    = True # DO_PLOT_GAIA
    DO_PLOT_NH_GAIA  = False
    DO_LABEL_NH_GAIA = False
    
    plot_tick_every  = 24*60*60
    hbt.figsize((15,10))
    mag_limit       = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    
    DO_PLOT_LORRI   = False
    DO_SCALEBAR_ARCSEC = True
    DO_UNITS_TITLE_DAYS = True
    DO_UNITS_TICKS_DAYS = True
    DO_PLOT_UNCERTAINTY_MU69 = True

#    DO_LABEL_USNO = False
#    DO_LABEL_GAIA = False
    
fov_lorri = 0.3*hbt.d2r  # Radians. LORRI width.

num_dt = 200 # Number of timesteps. 4000 is way too many -- 200 should be fine.

# Define plot colors

color_kbo        = 'purple'
color_kbo_radius = 'pink'
color_kbo_center = 'red'
color_lorri      = 'lightgreen'
color_roche      = 'yellow'
color_stars      = 'green'
color_gaia       = 'crimson'
color_tof        = 'blue'
color_uncertainty_mu69 = 'blue'

hbt.set_fontsize(fs)
    
try:
    lun = open(file_hd_pickle, 'rb')
    print("Loading file: " + file_hd_pickle)
    hd = pickle.load(lun)
    lun.close()
    
except IOError: # If Pickle file not found, then go ahead and read the catalog from scratch
    
    # Read HD into a table. Its positions are 1900. Coords are radians.

    print ("Reading HD catalog...")
    
    hd = read_hd_all()  
    num_stars = np.size(ra_1950)

    # Precess stars from 1900 to 2000
    
    print("Precessing to J2000")
    
    mx = sp.pxform('B1950', 'J2000', et[0]) # SPICE knows about 1950, but not 1900!
    
    ra_1950  = t['RA']
    dec_1950 = t['Dec']
        
    ra_2000  = []
    dec_2000 = []
    pt_1950  = []
    
    # Loop over every star and precess it individually. This loop is a real bottleneck.
    
    for i in range(num_stars):
        
        pt_1900 = sp.radrec(1, ra_1950[i], dec_1950[i])
        pt_1950 = sp.mxv(mx, pt_1900)
        pt_2000 = sp.mxv(mx, pt_1950) # SPICE will not precess 1900 -> 2000. So we apply 1950 -> 2000, do it 2x.
        d, ra_2000_i, dec_2000_i = sp.recrad(pt_2000)
        dec_2000.append(dec_2000_i)
        ra_2000.append(ra_2000_i)  
      
    hd['RA_2000']  = ra_2000
    hd['Dec_2000'] = dec_2000
    
    # Now save this pickle file, so we never have to run this routine again
    
    lun = open(file_hd_pickle, 'wb')
    pickle.dump(hd, lun)
    lun.close()
    print("Wrote: " + file_hd_pickle)

# Use the virtual observatory to get some UCAC2 stellar positions
# crval is a two-element array of [RA, Dec], in degrees
# Use conesearch.list_catalogs() to get list

encounter_phase = 'Inbound' if (et_start < et_ca) else 'Outbound'

# Set up the ET array

et = hbt.frange(et_start, et_end, num_dt)

index_asymptote = 0 if (encounter_phase == 'Outbound') else 1

# We make one call to SPICE to look up the asymptote position
# Most SPICE calls are later, but we need to do this one now, ahead of order.

(st,lt) = sp.spkezr('MU69', et[index_asymptote], 'J2000', 'LT', 'New Horizons')
(junk, ra_asymptote, dec_asymptote) = sp.recrad(st[0:3])   
crval = np.array([ra_asymptote, dec_asymptote]) * hbt.r2d
                  
radius_search = 0.15 # degrees # We only use these catalog for very fine searches, so narrow is OK.

DO_PLOT_GSC = False # Goes to about v=12

if DO_PLOT_GSC:                 
    name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
    stars = conesearch.conesearch(crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*hbt.d2r # Convert to radians
    dec_stars = np.array(stars.array['DEJ2000'])*hbt.d2r # Convert to radians
    mag_stars = np.array(stars.array['Pmag'])

if DO_PLOT_USNO:                 
    name_cat = u'The USNO-A2.0 Catalogue (Monet+ 1998) 1'
    stars = conesearch.conesearch(crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*hbt.d2r # Convert to radians
    id_stars = np.array(stars.array['USNO-A2.0']) # ID
    id_stars = id_stars.astype('U') # Convert from byte string to Unicode, ugh.
    dec_stars = np.array(stars.array['DEJ2000'])*hbt.d2r # Convert to radians
    mag_b_stars = np.array(stars.array['Bmag'])
    mag_r_stars = np.array(stars.array['Rmag'])

    usno = Table([id_stars, ra_stars, dec_stars, mag_b_stars, mag_r_stars], 
                 names = ['ID', 'RA_2000', 'Dec_2000', 'Bmag', 'Rmag'])

if DO_PLOT_GAIA:
    stars = gaia_query(crval[0], crval[1], radius_search, maxmag=20, maxsources=10000)
    ra_stars = np.array(stars['RA_ICRS'])*hbt.d2r # Convert to radians
    id_stars = np.array(stars['Source'])
    dec_stars = np.array(stars['DE_ICRS'])*hbt.d2r # Convert to radians
    mag_stars = np.array(stars['__Gmag_'])

    gaia = Table([id_stars, ra_stars, dec_stars, mag_stars], 
                 names = ['ID', 'RA_2000', 'Dec_2000', 'mag'])    

if DO_PLOT_NH_GAIA:
    stars = nh_gaia_query(crval[0], crval[1], radius_search, maxmag=20, maxsources=10000)
    ra_stars = np.array(stars['RA'])*hbt.d2r # Convert to radians
#    id_stars = np.array(stars['Source'])

    dec_stars = np.array(stars['Dec'])*hbt.d2r # Convert to radians
    
    # Do a bit of maniplation to take G mag, or B mag, or mean, depending on what we have
    
    mag_stars = np.array(stars['g'] + stars['r'])  
    is_both = (gaia['r'] * gaia['g']) > 0.0
    mag_stars[is_both] /= 2
             
    # Remove a very small number of stars that have clearly spurious magnitudes
         
    is_error = (mag_stars > 50)
    mag_stars[is_error] = 20          

    nh_gaia = Table([ra_stars, dec_stars, mag_stars], 
                 names = ['RA', 'Dec', 'mag'])  
    
#==============================================================================
# Look up NH position (ie, KBO position)
#==============================================================================

print("Looking up NH position...")


dt = et[1] - et[0] # Timestep

ra_kbo   = []
dec_kbo  = []
dist_kbo = []

ra_kbo_tof_minus = []
ra_kbo_tof_plus  = []
dec_kbo_tof_minus = []
dec_kbo_tof_plus  = []

for et_i in et:

    # Nominal position
    (st,lt) = sp.spkezr('MU69', et_i, 'J2000', 'LT', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3])   
    ra_kbo.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo.append(dec_i)
    dist_kbo.append(dist_i)

    # TOF uncertainty -- negative 
    (st,lt) = sp.spkezr('MU69', et_i-dt_tof, 'J2000', 'LT+S', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3]) 
    ra_kbo_tof_minus.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo_tof_minus.append(dec_i) 

    # TOF uncertainty -- positive 
    (st,lt) = sp.spkezr('MU69', et_i+dt_tof, 'J2000', 'LT+S', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3]) 
    ra_kbo_tof_plus.append(ra_i)
    dec_kbo_tof_plus.append(dec_i) 
  
# Convert all of these from lists to NP arrays
    
dec_kbo           = np.array(dec_kbo)
ra_kbo            = np.array(ra_kbo)
dec_kbo_tof_plus  = np.array(dec_kbo_tof_plus)
dec_kbo_tof_minus = np.array(dec_kbo_tof_minus)
ra_kbo_tof_plus   = np.array(ra_kbo_tof_plus)
ra_kbo_tof_minus  = np.array(ra_kbo_tof_minus)

ra  = hd['RA_2000']*hbt.r2d           # Look up stellar coordinates - be sure to use J2000 
dec = hd['Dec_2000']*hbt.r2d
dist_kbo = np.array(dist_kbo)*u.km  # NH-MU69 distance

# Compute the size of the Hill sphere

rho_kbo = 2.5*u.gram/u.cm**3
r_kbo   = 16.5 * u.km
m_kbo   = 4/3 * math.pi * r_kbo**3 * rho_kbo
a_kbo   = 43*u.au

a_hill = (a_kbo * (m_kbo / c.M_sun/3)**(1/3)).to('km')

# Compute angular size of Hill sphere
# In general the Hill sphere is just too large to plot... LORRI covers it only a couple days after C/A

angle_hill = (a_hill / dist_kbo).to('').value # Convert to radians, and drop the units
angle_hill = np.clip(angle_hill, -math.pi, math.pi) # Clip it (same as IDL > or < operator)

# Compute angular size of KBO itself

angle_kbo  = (r_kbo / dist_kbo).to('').value # Radians

# Compute angular size of Roche limit

r_roche     = 2.5 * r_kbo
angle_roche = (r_roche / dist_kbo).to('').value # Radians

# Set the x and y limits for the plot

xlim = hbt.mm(ra_kbo*hbt.r2d)   # Plotting limits, in deg
ylim = hbt.mm(dec_kbo*hbt.r2d)

xlim = np.array(xlim) + pad_ra_deg * np.array([-1,1]) # Pad the plot by the specified amount
ylim = np.array(ylim) + pad_dec_deg * np.array([-1,1])


# Filter the stars. Set a flag for each one we want to plot, based on mag and position.

is_good =  np.logical_and(
                            hd['Ptg'] < mag_limit,
               
                            np.logical_and( 
                                np.logical_and(xlim[0] < ra, 
                                               ra < xlim[1]),     
                                np.logical_and(ylim[0] < dec,
                                               dec < ylim[1])  )    )

# Generate human-readable strings for the plot

if (DO_UNITS_TITLE_DAYS):
    t_start_relative_str = "K{:+.0f}d".format((et_start - et_ca)/day)
    t_end_relative_str   = "K{:+.0f}d".format((et_end   - et_ca)/day)
else:
    t_start_relative_str = "K{:+}h".format((et_start - et_ca)/hour)
    t_end_relative_str   = "K{:+}h".format((et_end   - et_ca)/hour)     
    
#==============================================================================
# Make the plot
#==============================================================================

# Draw the main trajectory

fig, ax = plt.subplots()

DO_PLOT_TRAJECTORY = True

if DO_PLOT_TRAJECTORY:
    ax.plot(ra_kbo*hbt.r2d, dec_kbo*hbt.r2d, color=color_kbo_center)

# Plot the time-of-flight uncertainty.
# Concl: TOF error does not change the position of any of these occultations. It just changes the time 
# at which we see them. This wasn't obvious to me, but these two curves lay over each other identically.
    
if DO_PLOT_TOF:
    ax.plot(ra_kbo_tof_minus*hbt.r2d, dec_kbo_tof_minus*hbt.r2d, color=color_tof, marker = '+', ms=5)
    ax.plot(ra_kbo_tof_plus*hbt.r2d,  dec_kbo_tof_plus*hbt.r2d,  color=color_tof, marker = '+', ms=10)
    
# Draw the LORRI FOVs -- 0.3 x 0.3 deg square

if DO_PLOT_LORRI:
    dec_kbo_plus_lorri  = dec_kbo + (fov_lorri/2) * np.sqrt(2)  # All radians
    dec_kbo_minus_lorri = dec_kbo - (fov_lorri/2) * np.sqrt(2)
    
    
    ax.fill_between(ra_kbo*hbt.r2d, (dec_kbo_plus_lorri)*hbt.r2d, (dec_kbo_minus_lorri)*hbt.r2d, 
                     color=color_lorri, alpha=0.5, label = 'LORRI FOV')

# Draw the Hill sphere

ax.plot(ra_kbo*hbt.r2d, (dec_kbo - angle_hill)*hbt.r2d, linestyle = '--')
ax.plot(ra_kbo*hbt.r2d, (dec_kbo + angle_hill)*hbt.r2d, linestyle = '--')

# Draw the Roche radius

ax.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_roche)*hbt.r2d, (dec_kbo-angle_roche)*hbt.r2d, 
                 color=color_roche, alpha=0.5, label = 'MU69 Roche radius = {} km'.format(r_roche.to('km').value))

# Draw the MU69 radius

ax.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_kbo)*hbt.r2d, (dec_kbo-angle_kbo)*hbt.r2d, 
                 color=color_kbo_radius, alpha=1, label = 'MU69 radius = {} km'.format(r_kbo.to('km').value))

# Draw '+' symbols every hour

for i,et_i in enumerate(et):
    
    if (np.mod(et_i - et[0], plot_tick_every)) < dt:
        if (i == 0): # Pass the info for this to plt.label()... but only for the first call!
            kwargs = {'label' : 'MU69 position, SPICE, tcm22'}
        else:
            kwargs = {}
        ax.plot(ra_kbo[i]*hbt.r2d, dec_kbo[i]*hbt.r2d, marker='+', 
                 markersize=20, color='black', **kwargs )


ax.text(ra_kbo[0]*hbt.r2d, dec_kbo[0]*hbt.r2d, t_start_relative_str + '                  ' , 
        fontsize=12, clip_on=True, horizontalalignment='center')  # Why does plt.text() kill my plot??

ax.text(ra_kbo[-1]*hbt.r2d, dec_kbo[-1]*hbt.r2d, '                ' + t_end_relative_str, 
        fontsize=12, clip_on=True, horizontalalignment='center')  # Why does plt.text() kill my plot??

# Plot the HD stars
if DO_PLOT_HD:
    
    ax.plot(hd['RA_2000'][is_good]*hbt.r2d, hd['Dec_2000'][is_good]*hbt.r2d, linestyle='none', marker='.', 
         label = 'HD, V < {}'.format(mag_limit), color=color_stars)

# Label the HD stars

    if (DO_LABEL_HD): 

        for i in range(np.size(is_good)):
            if (is_good[i]):
                    string_hd = hd['ID'][i]
            else:
                string_hd = ''
            
            ax.text(hd['RA_2000'][i]*hbt.r2d, hd['Dec_2000'][i]*hbt.r2d, 
                '  {}  {:.1f}  {}'.format(hd['Type'][i], hd['Ptg'][i], string_hd),
                fontsize = 8, clip_on = True)

# Plot the USNO stars
#plt.plot(usno['RA_2000'][is_good]*hbt.r2d, hd['Dec_2000'][is_good]*hbt.r2d, linestyle='none', marker='.', 
#         label = 'HD star', color=color_stars)

if DO_PLOT_USNO:
    ax.plot(usno['RA_2000']*hbt.r2d, usno['Dec_2000']*hbt.r2d, linestyle='none', marker='.', 
         label = 'USNO positions', color=color_stars)

if DO_LABEL_USNO:
    for i in range(np.size(usno)):
        ax.text(usno['RA_2000'][i]*hbt.r2d, usno['Dec_2000'][i]*hbt.r2d, 
            '  B={:.1f}, R={:.1f}  {}'.format(usno['Bmag'][i], usno['Rmag'][i], usno['ID'][i]),
            fontsize = 8, clip_on = True)

# Plot the Gaia stars

if DO_PLOT_GAIA:
    ax.plot(gaia['RA_2000']*hbt.r2d, gaia['Dec_2000']*hbt.r2d, linestyle='none', marker='.', 
         label = 'Gaia positions', color=color_gaia, alpha = 0.5)

if DO_LABEL_GAIA:
    for i in range(np.size(gaia)):
        ax.text(gaia['RA_2000'][i]*hbt.r2d, gaia['Dec_2000'][i]*hbt.r2d, 
            '  v={:.1f}  {}'.format(gaia['mag'][i], gaia['ID'][i]),
            fontsize = 8, clip_on = True)
    
# Plot a scalebar if requested

if (DO_SCALEBAR_ARCSEC):

    d2as = 60. * 60.
    as2d = 1/d2as
    y0 = ylim[0]
    dy = ylim[1] - ylim[0]
    
    delta_pos_usno_as = 0.2
    delta_pos_mu69_as = 0.5 # This is positional uncertainty now, from Earth. I want to map this into 
    
    # Get the distance to Pluto for start and end times
    
    y0_scalebar = ylim[0] + 0.1*(ylim[1]-ylim[0])
    x0_scalebar = xlim[0] + 0.1*(xlim[1]-xlim[0])
    x1_scalebar = x0_scalebar + 1 * as2d
 
    ax.hlines(y0 + dy * 0.1, x0_scalebar, x0_scalebar + 1*as2d)
    ax.hlines(y0 + dy * 0.15, x0_scalebar, x0_scalebar + delta_pos_usno_as*as2d)
#    plt.hlines(y0 + dy * 0.2, x0_scalebar, x0_scalebar + delta_pos_mu69_as*as2d)
    
    ax.text(x0_scalebar, y0_scalebar,  "    = 1\" = {:.0f} km (at {}) = {:.0f} km (at {})".format(
                                                                          (dist_kbo[ 0].value)*1*hbt.as2r,
                                                                          t_start_relative_str,
                                                                          (dist_kbo[-1].value)*1*hbt.as2r,
                                                                          t_end_relative_str ), 
                                                         bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
#    ax.figtext(0.23, 0.23, "USNO accuracy = {}\"".format(delta_pos_usno_as))
#    plt.figtext(0.23, 0.265, "MU69 uncertainty = {}\"".format(delta_pos_mu69_as))

# Plot error ellipses for the position of MU69. These are from JS's email to me, which are in turn from Marc Buie's
# slides as of 2-Mar-2017. I am assuming that everything is equatorial.


if (DO_PLOT_UNCERTAINTY_MU69):
    
    dpos_x = 6413*u.km  # Uncertainty in X position, km, from JS email 22-Mar-17. Halfwidth
    dpos_y = 366*u.km   # Uncertainty in Y position, km

    angle = 3           # Rotation angle of ellipse, in degrees. This is a guess, just to make 
                         # it vaguely align with what is in plots of Buie (from JS email).
    
    # Plot at start of period
    
    width  = (dpos_x/dist_kbo[0]).value*hbt.r2d*2/math.cos(dec_kbo[0])
    height = (dpos_y/dist_kbo[0]).value*hbt.r2d*2
    
    xy = (ra_kbo[0]*hbt.r2d, dec_kbo[0]*hbt.r2d)  # Get uncertainty in x and y, and convert to deg
    ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, angle = angle, alpha=0.1, 
                                     color=color_uncertainty_mu69, label = 'MU69 3$\sigma$ pos, Buie')
    ax.add_patch(ell)    

    # Plot at end of period
    
    width  = (dpos_x/dist_kbo[-1]).value*hbt.r2d*2/math.cos(dec_kbo[-1])
    height = (dpos_y/dist_kbo[-1]).value*hbt.r2d*2
    
    xy = (ra_kbo[-1]*hbt.r2d, dec_kbo[-1]*hbt.r2d)  # Get uncertainty in x and y, and convert to deg
    ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, angle = angle, alpha=0.1, 
                                     color=color_uncertainty_mu69)
    ax.add_patch(ell)    


#==============================================================================
# Plot a circle showing where 1000 km radius rings would be, if they were there
#==============================================================================

DO_PLOT_RING_MU69 = True

if DO_PLOT_RING_MU69:
    
    dpos_x = 3000*u.km   # Radius of the rings to draw
    dpos_y = 3000*u.km
    
    angle = 0           # Rotation angle of ellipse, in degrees
    
    # Plot ring at start of period
    
    width = (dpos_y/dist_kbo[0]).value*hbt.r2d*2/math.cos(dec_kbo[0])
    height = (dpos_y/dist_kbo[0]).value*hbt.r2d*2
    
    xy = (ra_kbo[0]*hbt.r2d, dec_kbo[0]*hbt.r2d)  # Get uncertainty in x and y, and convert to deg
    ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, angle = angle, alpha=0.5, 
                                     facecolor='none', edgecolor = 'grey', linewidth=3)
    ax.add_patch(ell)    

    # Plot ring at end of period
    
    width = (dpos_y/dist_kbo[-1]).value*hbt.r2d*2/math.cos(dec_kbo[-1])
    height = (dpos_y/dist_kbo[-1]).value*hbt.r2d*2
    
    xy = (ra_kbo[-1]*hbt.r2d, dec_kbo[-1]*hbt.r2d)  # Get uncertainty in x and y, and convert to deg
    ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, angle = angle, alpha=0.5, 
                                     edgecolor='grey', facecolor='none', linewidth=3, 
                                     label = 'MU69 ring, radius = {:.0f} km'.format(dpos_x.to('km').value))
    ax.add_patch(ell)

#==============================================================================
# Plot a Lauer image with MU69 at center, and axes in km.
#==============================================================================
#%%

radius_ring = 3000*u.km
radius_plot = 7400  # Radius, in km

color_stars = 'blue'

fig, ax = plt.subplots()

# Plot ring around MU69


ax.set_xlim(radius_plot * np.array([-1,1]))
ax.set_ylim(radius_plot * np.array([-1,1]))

# Grab RA and Dec for all stars

ra_stars  = gaia['RA_2000'] 
dec_stars = gaia['Dec_2000']

for i,t in enumerate(et[0:1000]):  # Plot for every timestep

    d_ra_ang  = (ra_stars  - ra_kbo[i])/np.cos(dec_stars)
    d_dec_ang = dec_stars - dec_kbo[i]
  
    d_ra_km_proj = d_ra_ang * dist_kbo[i]
    d_dec_km_proj = d_dec_ang * dist_kbo[i]

    ax.plot(d_ra_km_proj, d_dec_km_proj, marker = '.', linestyle='none', color = color_stars, markersize=1)

    # Label each star at its starting position (at outer edge)

    if (i == 0):
        for j in range(np.size(gaia)):  # Loop over stars
            ax.text(d_ra_km_proj[j].value, d_dec_km_proj[j].value, '  mag={:.1f}'.format(gaia['mag'][j]),
                fontsize = 8, clip_on = True)
        ax.plot(d_ra_km_proj[i], d_dec_km_proj[i], marker = '.', linestyle = 'none', color=color_stars,
                markersize=3, label = 'Gaia star')
#            print("Plotted at {}, {}, {}".format(d_ra_km_proj[j], d_dec_km_proj[j], gaia['mag'][j]))
        
# Plot MU69

ax.plot(0,0,marker = 'o', color = color_kbo, label = 'MU69', linestyle='none')
ax.set_xlabel('Projected Distance, RA [km]')
ax.set_ylabel('Projected Distance, Dec [km]')

ax.set_title('{} .. {}'.format(t_start_relative_str, t_end_relative_str))
    
xy = np.array([0,0])
width  = radius_ring.to('km').value * 2
height = radius_ring.to('km').value * 2
angle  = 0

ell = matplotlib.patches.Ellipse(xy = xy, width=width, height=height, angle = angle, alpha=0.5, 
                                     edgecolor='grey', facecolor='none', linewidth=3, 
                                     label = 'MU69 ring, radius = {:.0f} km'.format(radius_ring.to('km').value))
ax.add_patch(ell)
ax.set_aspect('equal')

ax.legend(loc = 'lower right')    
plt.show()
plt.savefig('up.png')

# Calculate the LORRI pixel scale at start and end of this.

ang_pix_lorri = 0.3*hbt.d2r/1024 # Radians per pixel, LORRI

km_pix_lorri_start = ang_pix_lorri * dist_kbo[0].to('km').value
km_pix_lorri_end   = ang_pix_lorri * dist_kbo[-1].to('km').value

print("At {}, LORRI pixel scale = {:.0f} km/pix; plot width = {:.0f} pix".format(
        t_start_relative_str, km_pix_lorri_start, 2*radius_plot/km_pix_lorri_start))
print("At {}, LORRI pixel scale = {:.0f} km/pix; plot width = {:.0f} pix".format(
        t_end_relative_str, km_pix_lorri_end, 2*radius_plot/km_pix_lorri_end))
print("")

#print("Plot width = {:.0f} LORRI pixels start, {:.0f} pixels end".format(
#        2*radius_plot/km_pix_lorri_start, 2*radius_plot/km_pix_lorri_end))

#%%
    
# Convert RA and Dec to to projected distance at this et  
# Assume a ring size 
# Plot every time that a star goes within, say, 100 km 
    
#fig, ax = plt.subplots()
    
#mean = [ 0.5 ,  0.5]
#width = 0.1
#height = 0.2
#angle = -54
#fig, ax = plt.subplots()
#plt.plot(x,y)

    
       
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')

  
        
if (DO_UNITS_TITLE_DAYS):
    ax.set_title('MU69 {}, ticks every {:.0f}h, {} .. {}'.format(
              encounter_phase, plot_tick_every/(60*60), t_start_relative_str, t_end_relative_str))
else:    
    ax.set_title('MU69 {}, ticks every {:.0f} min, K{} .. {}'.format(
              encounter_phase, plot_tick_every/60, t_start_relative_str, t_end_relative_str))
    
ax.legend(framealpha=0.8, loc = 'lower right')
plt.show()

#==============================================================================
# Do some test cases
#==============================================================================

# Print coords for one candidate star

#for i in range(np.size(usno)):             # Ugh -- loop over and check each one!
#                                           # Should be able to use np.equal() but 
#                                           # diagnostic says it is not actually implemented yet?!
#    if (usno['ID'][i] == '1050-03305028'): # This is the good occultation star
#        print("ID={}, RA={} deg, Dec={} deg".format(
#          usno['ID'][i],
#          usno['RA_2000'][i]*hbt.r2d,
#          usno['Dec_2000'][i]*hbt.r2d))
#      
 
#==============================================================================
# Make a plot of USNO vs. Gaia stellar positions. Run this after the variables are initialized in above code.
#==============================================================================

DO_TEST_GAIA = False

if DO_TEST_GAIA:
    plt.plot(gaia['RA_2000']*hbt.r2d, gaia['Dec_2000']*hbt.r2d, marker = '.', linestyle='none', color = 'blue', 
             label = 'Gaia positions')
   
    plt.plot(usno['RA_2000']*hbt.r2d, usno['Dec_2000']*hbt.r2d, linestyle='none', marker='o', 
         label = 'USNO positions', color='none', markeredgecolor='green', mew=1)
    plt.plot()
    plt.xlabel('RA [deg]')
    plt.title('Gaia vs USNO positions')
    plt.ylabel('Dec [deg]')
    plt.plot(crval[0], crval[1], marker = 'o', color = 'pink', ms=20, alpha=0.5, linestyle='none', 
             label='MU69 asymptote')
    DO_WINDOW_ZOOM = True
    if (DO_WINDOW_ZOOM):
        plt.xlim((274.70, 274.80))
        plt.ylim((-20.95, -20.80))
    plt.legend(loc='upper right')
    plt.show()
    
    