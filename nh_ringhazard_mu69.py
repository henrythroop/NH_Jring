"""
     Calculate the hazard to the NH spacecraft based on the I/F or tau
     This code generates the plots that go into my NH Ring Hazard whitepaper, Nov 2011.
     
     Adapted from ~/occ/AAT/HBT/nh_ringhazard_v2_bw.pro
     
      HBT Nov 2011 Original version
"""      
                        
#; HBT Nov 2012 Updated to include phase function explicitly in backscatter calculation as per MMH.
#;   In earlier version, I was explicitly assuming isotropic scatterers because I did not include
#;   phase function P11. This was a mistake: for big grains, most of the photons end up into the 
#;   diffraction peak, and the q_sca asymptotes to 0.5 to show this. However, the fraction of photons
#;   getting back to the observer is essentially zero, because such large Mie particles are very
#;   absorptive. The value of q_sca * p11 approaches zero for backscatter, and that is what I should
#;   have used instead. As a result, the visibility of big grains was under-estimated (by 10x or more),
#;   and the number N of dangerous grains was similar under-estimated.
#;
#;   However, this whole discussion may be academic: 200 um grains (X=2500) are huge, and probably 
#;   Mie scattering is a terrible approximation for them. If they are cracked, rough, non-spherical,
#;   stuck to one-another, heterogeneous, fractured... all of these will cause deviations from Mie, 
#;   which will crank up the effective P11 at backscatter, making big grains more visible, and thus
#;   lowering N back down. My original white paper didn't discuss non-spherical grains, because it's 
#;   very hard to get any handle on what these things might look like. 
#
#;   Mie scattering says that hailstones should be invisible (that is, near-perfect forward-scatterers),
#;   but experience shows otherwise!
#; 
#;   NB: In the output files, 'v1' indicates original version of these files, kept around for historical
#;       reasons. 'v2' means the revised and fixed files.
#;  *BW* is a new version of this code, identical except it outputs BW files instead of color
		

# Q: how to best handle the upper limit of the size distribution?
# o Cut it off at the size at which there will be *one grain* of that size (like 10 meters)
# o Look at the Gruen GZF84 size dist, and pick a flux to cut off at, and then cut our size off there basically.
#   Looks like this woiuld be a cutoff at something like 1d-4 .. 1d-0 grams. Let's say the largest impactor is a gram.
#   And the yield is 1d-8 (m_ej = 1d-8 E_imp, in cgs), P. 11 GWH78. 

# o Cut it off arbitrarily at, say, 
# o Cut it off at the size where the lifetime (against radiation pressure, drag, etc) is the age of the ring system. 
#   I could get an estimtae for this from Hamilton.

# Changes to make for MU69:
# Put in terms of I/F. That is what most of our     


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

from   astropy.utils import data
from   astropy.wcs import WCS

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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
#import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

from   matplotlib.figure import Figure

from   pymiecoated import Mie

# HBT imports

import hbt

# Pick a particle scattering model. We can use Mie, or Lambertian

DO_MIE      = True
DO_LAMBERT  = True
do_merged   = True			# Mie and Lambert, mixed and cross-faded

DO_HST = False
DO_AAT = False

DO_MU69_INBOUND = True

pi     = math.pi

# Write n(r) out to a file? For Hal Weaver

do_write_txt_file = True

# Set the wavelength 'alam'

alam		= 500 * u.nm
  
# Set the plot thickness

linethick     = 3

# Set the transition radius between the two scattering models. Same as in TE98.

r_trans	= 300 * u.nm
  
halfwidth_trans	= 3		# transition goes from fraction below r_trans, to this multiple above it

# Set the phase angle of the HST observations

phase_inbound = 10.5*u.deg   # Phase angle, NH, inbound MU69  (from GV)
phase_outbound = 168.2*u.deg # Phase angle, NH, outbound MU69 (from GV)
phase_hst      = 2.0 * u.deg
 
# Define the albedo. This was not used in v1. It is used in v2, as a limiting case.

albedo	= 0.05			# We pick the best-guess, which is Charon's albedo = 0.35. 0.05 is more conservative,
  					# but may be meaninglessly low.

# Set the spacecraft area

sarea = 5*(u.m)**2     # SS07 used 10 m^2. But SAS e-mail 27-Nov-2011 says to use 5.5 m^2 instead. Hersman sayd 3.5 m^2. 
                    # So I call it 5.

# Set the particle size range

r_danger_1 = 0.2*u.mm   # This is the size defined as the critical size by SAS in Nov-2011 meeting
  
nhits_1 = 0.1
nhits_2 = 0.01

rmin   = 0.01 * u.micron
rmax   = 5000 * u.micron   # 900 microns is around the cutoff size of the UK Mie code (X ~ 12,000).
                    # By coincidence, this is also very close to the maximum particle size that 
                    # I calculate for the ring (see code attached, below)

num_r    = 50       # Number of radial bins in Mie calculation

iof_norm = 2e-10
tau_norm = iof_norm
phase = phase_inbound
ylim = (1e-15, 1e6)
title = r'MU69 Inbound'
    
alpha        = pi*u.rad - phase	    # Alpha = scattering angle, in degrees
cos_alpha    = np.cos(alpha)
dqv          = cos_alpha.value              # .value -- take it out of 'units' world
    
subobslat_nh  = -49.07*u.deg                 # 2015 Mar 26, Pluto from NH

# Create the plot

q         = np.array([1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7])			# 9 elements long
colors    = np.array(['black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black'])
linestyle = np.array(['-', ':', '--', '-.', '-:', '__', '-', ':', '--'])
linethick = np.array([1, 1, 1, 1, 1, 1, 3, 3, 3])
  
#p = plot([3,4], xtitle = 'Radius [$\mu$m]' ,$
#      ytitle = 'Total number N > r impacting NH during encounter',xrange=[0.1,1000d],yrange=yrange, 
# /XLOG, /YLOG, /NODATA, $
#      title = title)
   
#p_arr = replicate(p, sizex(q)) 

# Draw a dotted line at the two danger sizes that SAS specified at 3-Nov-2011 workshop

r    =  hbt.frange(rmin, rmax, num_r, log=True)*u.micron  # Astropy bug? When I run frange(), it drops the units.
qmie = np.zeros(num_r)
qsca = np.zeros(num_r)
qext = np.zeros(num_r)
qbak = np.zeros(num_r)
qabs = np.zeros(num_r)
p11_mie  = np.zeros(num_r)

bin_r_danger_1 = np.where(r > r_danger_1)[0][0]

# Calculate the transition regime

binstart_trans	= hbt.wheremin(np.abs(r - r_trans/halfwidth_trans))
  
binend_trans	    = hbt.wheremin(np.abs(r - r_trans*halfwidth_trans))

frac_mie	= np.concatenate([1 + np.zeros(binstart_trans), 
                         hbt.frange(1, 0, binend_trans-binstart_trans), 
                         0 + np.zeros(num_r-binend_trans)])
frac_lambert	= 1 - frac_mie

# Calc Q_ext *or* Q_sca, based on whether it's reflected or transmitted light
 
n_refract = 1.33
m_refract = -0.001
nm_refract = complex(n_refract,m_refract)
x = (2*pi * r/alam).to('1').value
   
print('Doing Mie code')

# Mie code doesn't compute the phase function unless we ask it to, by passing dqv.
 
if DO_MIE:
    for i,x_i in enumerate(x):  # Loop over particle size X, and get Q and P_11 for each size
        mie = Mie(x=x_i, m=nm_refract)  # This is only for one x value, not an ensemble
        qext[i] = mie.qext()
        qsca[i] = mie.qsca()
        qbak[i] = mie.qb()  
      
        (S1, S2)  = mie.S12(np.cos(pi*u.rad - phase)) # Looking at code, S12 returns tuple (S1, S2).
                                                             # For a sphere, S3 and S4 are zero.
                                                             # Argument to S12() is scattering angle theta, not phase
        k = 2*pi / alam
      
        sigma = pi * r[i]**2 * qsca[i] # Just a guess here
      
# Now convert from S1 and S2, to P11: Use p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf
      
        p11_mie[i]  = 4 * pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2
       
# Now assume a I/F = 1. For the current size dist, if I/F = 1, then what # of impacts will we have?

              
p11_lambert		 = 8/(3 * pi) * (np.sin(phase) + (pi - np.sin(phase)) * np.cos(phase)).value	
  
  			          # But hang on -- how is this phase function right since it doesn't include diffraction! Ugh...
  			          # The real phase function should really be half this, I think. This properly integrates to 
  			          # 2 over 0 .. 2 pi -- but what *should* integrate to 2 is the phase function including 
  			          # the diffraction spike, meaning these P11's in backscatter would go down.
  			          # OK, halving it. That will have the effect of making particles darker, harder to see, and 
  			          # increasing N for lambert.
  
p11_lambert      *= 0.5
  			          
p11_lambert      *= albedo

qsca_mie          = qsca
qsca_lambert		= 2	# I think this is right. Scattering efficiency = 2 since it can scatter and diffract?

p11_merged		= (p11_mie  * frac_mie) + (p11_lambert  * frac_lambert)
qsca_merged		= (qsca_mie * frac_mie) + (qsca_lambert * frac_lambert)

print('Done UK Mie code')

if (DO_AAT):
    qmie = qext
      
if (DO_HST):
    qmie = qsca

# Draw a polygon for the filled dangerous area

# Set a large font size

matplotlib.rc('font', size=15)

# Set up the two subplots

hbt.figsize((10,14))
(fig, axes) = plt.subplots(2, 1)

# Loop over each value of q

for q_i in q:
  
    q_nominal = 2.5
    
    if (q_i == q_nominal):
        linewidth = 4
    else:
        linewidth = 2
        
    n    = hbt.powerdist(r, 1, q_i, 1).value
  
# Now normalize this to be number per cm2, in cross-section.
 
    n_lambert = (n * tau_norm / np.sum(n * pi * r**2 * qsca_lambert * p11_lambert)).to('1/cm^2')    # Added P11 20-Nov-2012
    n_mie     = (n * tau_norm / np.sum(n * pi * r**2 * qsca_mie     * p11_mie)).to('1/cm^2')        # Added P11 20-Nov-2012
    n_merged  = (n * tau_norm / np.sum(n * pi * r**2 * qsca_merged  * p11_merged)).to('1/cm^2')

# Now convert into total hitting s/c -- that is, N > r, for many different values of r.
# Also, take out of units system

    n_cum_sc_lambert = (hbt.incr2cum(n_lambert, DO_REVERSE=True) * sarea / np.cos(subobslat_nh)).to('1').value
    n_cum_sc_mie     = (hbt.incr2cum(n_mie,     DO_REVERSE=True) * sarea / np.cos(subobslat_nh)).to('1').value
    n_cum_sc_merged  = (hbt.incr2cum(n_merged,  DO_REVERSE=True) * sarea / np.cos(subobslat_nh)).to('1').value

# Now assume I/F = 1. Given this, how many particles would we hit (of any size), in surface area of s/c?

    iof_reference = 1 # Define the reference I/F
    
    iof_default = np.sum(n_merged * pi * r**2 * qsca_merged * p11_merged) * albedo / 4
    
                        # This comes from white paper eq1, or TPW eq1.
                        
# Calculate the n (# cm^-2) necessary to achieve iof = 1
                        
    n_for_iof_reference = n_merged / iof_reference  # Number cm^-2 of total particles, for I/F = 1
    
    num_iof = 100 #  Number of I/F bins to use
    
    iof_arr = hbt.frange(1e-10, 1, num_iof, log=True)

    n_damaging_for_iof_reference = np.sum(n_for_iof_reference[bin_r_danger_1:]) #  Number cm^-2 damaging, for I/F = 1
    N_damaging_for_iof_reference = (n_damaging_for_iof_reference * sarea).to('1').value
                                               #  Total number damaging, that hit sc, for I/F = 1
    
    N_damaging_arr = (iof_arr / 1) * N_damaging_for_iof_reference # Total num damaging sc, for all I/F

# XXX Redo calculation above from scratch. I want to get # of fatal hits, for this distribution

    n_fatal = np.sum(n_merged[bin_r_danger_1:])
    
# Add a line for (radius) vs (# of hits), for this q, on Plot #1

    axes[0].plot(r.to('micron').value, n_cum_sc_merged, label = 'q = ' + repr(q_i), linewidth=linewidth)

# Add a line for (I/F) vs (# of hits), for this q, on Plot #2

    axes[1].plot(iof_arr, N_damaging_arr, label = 'q = {}'.format(q_i), linewidth=linewidth)
    axes[1].set_yscale('log')
    axes[1].set_xscale('log')
    axes[1].set_ylabel(r'N > $r_c$ during passage')
    axes[1].set_xlabel(r'Observed I/F limit at $\lambda$ = {} $\mu$m'.format(alam.to('micron').value))
    axes[1].legend()
    
# Label r_c

axes[0].axvline(r_danger_1.to('micron').value,linestyle=':')

# Draw dotted lines at the number of hits that SAS specified

axes[0].axhline(nhits_1, linestyle = ':')
axes[0].axhline(nhits_2, linestyle = ':')

axes[1].axhline(nhits_1, linestyle=':')
axes[1].axhline(nhits_2, linestyle=':')

# Done with plotting the individual lines... now set the axis limits, and finish and render the plots

axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlabel('Radius [$\mu$m]')
axes[0].set_ylabel(r'Total N > $r_c$ impacting NH during encounter')
axes[0].set_xlim((0.1, 1000))
axes[0].set_ylim(ylim)
axes[0].set_title(title + ", I/F = {}".format(iof_norm))
axes[0].legend()


axes[1].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_ylim((1e-24, 10))
axes[1].set_xlabel('Observed dust I/F limit')
axes[1].set_ylabel(r'Total N > $r_c$ impacting NH during encounter')
axes[1].set_title(title + r', $r_c$ = {:.0f} $\mu$m, a = {}'.format(r_danger_1.to('micron').value, albedo))
axes[1].legend()

fig.show()

# Now make a plot showing (# of hits) vs (I/F) for each of these curves

# Draw the word 'Danger' in the UR quadrant

#  te = text(r_danger_1 * 1.3/um, ((p.yrange)[1])/100, 'Danger', /DATA, /OVERPLOT)

# Draw the legend

#  pos_legend = [1,5d-3]        ; Upper-left corner of legend is put at this location
  
     
#==============================================================================
#  Verify normalization of phase functions P11, etc.
#==============================================================================

def other():

    theta_0 = 0
    theta_1 = pi
    num_theta = 100
    
    theta = hbt.frange(theta_0, theta_1, num_theta)
    dtheta = theta[1] - theta[0]
    
    p11 = np.cos(theta)  		# lambertian phase function = cos(theta). This is the *surface* phase function --
      				# not the phase function for the entire sphere, which needs to be convolved with the
    				# crescent shape!
    
    p11 = 1			# Integral = 2 for this value (isotropic scattering)
    p11 = 8/(3 * pi) * (np.sin(theta) + (pi - theta) * np.cos(theta))			
      				# Integral = 2 for this expression (disk-integrated scattering from lambert sphere)
      				# Madhusudhan & Burrows, 2012, http://arxiv.org/pdf/1112.4476v1.pdf, eq33
    
    int = np.sum(p11 * np.sin(theta) * dtheta)		# This integral should be 2, as per TE98 eq 16.
 
    print("Total = {}".format(int))

#==============================================================================
# Calculate the largest particle size in the ring.
#==============================================================================

# The philosophy here is that according to Gruen et al 1985, the mass flux of interplanetary dust takes a steep dropoff
# above 1 gram (cf fig. 3, GZF85, p. 10). So, I assume that a 1 gram particle is the largest that will hit  Pluto ring
# system. Maybe this is a bit ad hoc... their figure three shows a knee, but that may be just due to how it's plotted, 
# the real power dist may continue to go on forever.
# Flux of 1 gram particles is 2.2d-15/m2/sec at 1 AU. Maybe it's down by 100x at Pluto (my guess).
# That's basically one per day hitting Charon.

#  v_imp = 5 * km
#  m_imp	= 1d * gram
#  k_ej 	= 1d-8
#  m_ej 	= 0.5 * m_imp * v_imp^2 * k_ej
#  r0	= 0.01 * um
#  r1 	= 100*cm
#  q_ej  = 3.5
#  num_r	= 100
#  r	= frange(r0, r1, num_r, /log)
#  n	= powerdist(r, m_ej, q_ej, rho)	; Create the ejecta. Yes, most of this will not stay in the ring, but we are
#  						; conservative and ignore that.
#  rho	= 1
#  m	= total(4/3d * pi * r^3 * rho)
#
#  n	= n * (m_ej / m)			; Normalize the ejecta number to the proper mass
#  r_max	= r[max(where(n gt 1))]
#
#  print, 'For impactor size ' + st(m_imp) + ' grams, the max ejecta size is ' + st(r_max/um) + ' um.'

                                   
def test_p11():
    '''
    Test to plot a phase function, and make sure it is normalized properly
    '''
    
    alam = 500*u.nm
    r = 10*u.micron
    x = (2*math.pi * r / alam).to('1').value
    

    num_theta = 1000
    theta = hbt.frange(0,math.pi, num_theta)
    p11 = np.zeros(num_theta)
    
    phase = math.pi - theta      # Theta is scattering angle
    
    nm_refract = complex(1.5, 0.1)
    mie = Mie(x=x, m=nm_refract)  # This is only for one x value, not an ensemble

    qext = mie.qext()
    qsca = mie.qsca()
    qbak = mie.qb()

    for i,theta_i in enumerate(theta):
        (S1, S2)  = mie.S12(np.cos(theta_i)) # Looking at code, S12 returns tuple (S1, S2). S3, S4 are zero for sphere.
        k = 2*pi / alam
        sigma = pi * r**2 * qsca             # For (S1,S2) -> P11: p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf
                                             # qsca = scattering efficiency. sigma = scattering cross-section 
        p11[i]  = 4 * pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2
                   
    # Check the normalization of the resulting phase function.

    dtheta = theta[1] - theta[0]
    
    norm = np.sum(p11 * np.sin(theta) * dtheta)  # This should be 2, as per TPW04 eq. 4
    
    print('Normalized integral = {:.2f}'.format(norm))
    
    plt.plot(phase*hbt.r2d, p11)
    plt.yscale('log')
    plt.title('X = {:.1f}'.format(x))
    plt.xlabel('Phase angle')
    plt.ylabel('$P_{11}$')
    plt.show()