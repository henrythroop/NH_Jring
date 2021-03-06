"""
     Calculate the hazard to the NH spacecraft based on the I/F or tau
     This code generates the plots that go into my NH Ring Hazard whitepaper, Nov 2011.
     
     Adapted from ~/occ/AAT/HBT/nh_ringhazard_v2_bw.pro
     
     ** This file will create exactly the same plot as in nh_ringhazard_v2.pro (the IDL version). 
        It is held here as a baseline version, just to compare IDL vs. Python.
     
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
from   astroquery.vo_conesearch import conesearch
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
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

DO_HST = True
DO_AAT = not(DO_HST)

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

phase_hst	= 2.0*u.deg	# Phase angle of HST observations, in degrees.
  					# GV calculates range of 0 .. 2 deg for Earth-Pluto phase angle in summer 2011.
					# Actual results are insensitive to the precise angle because phase function of both 
					# small and large grains are flat in this region.
					# For HST observations, Phase ~ 1 deg -> scattering ~ 179 deg.
					# MIE_SINGLE requires the scattering angle, not the phase angle.

alpha_hst	= pi*u.rad - phase_hst	    # Alpha = scattering angle, in degrees
cos_alpha 	= np.cos(alpha_hst)
dqv		     = cos_alpha.value                 # .value -- take it out of 'units' world

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
num_r    = 500

if (DO_AAT):
    tau_norm  = 0.0035    # This is the input optical depth. Determined from the occultation, from AAT_hbt.pro (part 7).
    ylim = (1e-10, 1e10)
    file_out = 'pluto_dust_limits_aat_v2_rmax' + repr(rmax.to('micron')) + '_bw'
    title = r'Pluto dust limits, 2006 AAT occultation, $\tau$ < 0.0035'  # Use 'r' to denote raw string, not escaped

if (DO_HST):
    iof_norm	= 3e-7
    tau_norm	= iof_norm
    ylim = (1e-15, 1e6)
    file_out = 'pluto_dust_limits_hst_v2_rmax' + repr(rmax.to('micron')) + '_bw'
    title = r'Pluto inner region, HST 2012 imaging, I/F < 3$\times 10^{-7}$, $r_{max}$ = ' + \
        repr(rmax.to('micron').value) + ' $\mu$m'

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
    for i,x_i in enumerate(x):  
        mie = Mie(x=x_i, m=nm_refract)  # This is only for one x value, not an ensemble
        qext[i] = mie.qext()
        qsca[i] = mie.qsca()
        qbak[i] = mie.qb()  
      
        (S1, S2)  = mie.S12(np.cos(pi*u.rad - phase_hst)) # Looking at code, S12 returns tuple (S1, S2).
                                                             # For a sphere, S3 and S4 are zero.
                                                             # Argument to S12() is scattering angle theta, not phase
        k = 2*pi / alam
      
        sigma = pi * r[i]**2 * qsca[i] # Just a guess here
      
      # Now convert from S1 and S2, to P11: Use p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf


#        (S1, S2)  = mie.S12(np.cos(theta_i)) # Looking at code, S12 returns tuple (S1, S2). S3, S4 are zero for sphere.
#        k = 2*pi / alam
#        sigma = pi * r**2 * qsca             # For (S1,S2) -> P11: p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf
                                             # qsca = scattering efficiency. sigma = scattering cross-section 
#        p11[i]  = 4 * pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2


      
        p11_mie[i]  = 4 * pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2
      
#      mie_single, x, nm_refract, dqxt, dqsc, dqbk, dg, xs1, xs2, dph, dqv=dqv 
  							# UK Mie code. m_refract must be neg. x must be < 12,000.
                              # XXX I don't get this. Is X a scalar, or vector?
       
#    qext = dqxt		# Q_ext
#  qsca = dqsc		# Q_sca
#  qbak = dqbk		# Q_back (??)
#  p11_mie  = dph[*]	# Phase function. Computed only if angles dqv is passed.
 
p11_lambert		 = 8/(3 * pi) * (np.sin(phase_hst) + (pi - phase_hst) * np.cos(phase_hst)).value	
  
                        # Nov-2017. Corrected lambertian calculation. Was using "(pi-np.sin(phase_hst))".
                        # Shouldn't make any difference since we are at backscatter and second term is ~0.
                        
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

#  po = polygon([r_danger_1/um,(p.xrange)[1],(p.xrange)[1],r_danger_1/um], $
#    [nhits_1, nhits_1,(p.yrange)[1],(p.yrange)[1]],/data, $
#    fill_color='pink',fill_transparency=0,/over, linestyle=' ')

for q_i in q:		# Loop downward, so that the low q's are plotted on top for visibility
  
    n    = hbt.powerdist(r, q_i, mass=1, rho=1).value
  
# Now normalize this to be number per cm2, in cross-section.
 
    n_hst_lambert = n * tau_norm / np.sum(n * pi * r**2 * qsca_lambert * p11_lambert) 	# Added P11 20-Nov-2012
    n_hst_mie     = n * tau_norm / np.sum(n * pi * r**2 * qsca_mie     * p11_mie) 	# Added P11 20-Nov-2012
    n_aat         = n * tau_norm / np.sum(n * pi * r**2 * qsca_mie)
    n_hst_merged  = n * tau_norm / np.sum(n * pi * r**2 * qsca_merged * p11_merged)

# Now convert into total hitting s/c
# Also, take out of units system

    n_cum_sc_hst_lambert = (hbt.incr2cum(n_hst_lambert, DO_REVERSE=True) * sarea / np.cos(subobslat_nh)).to('1').value
    n_cum_sc_hst_mie 	   = (hbt.incr2cum(n_hst_mie, DO_REVERSE=True)     * sarea / np.cos(subobslat_nh)).to('1').value
    n_cum_sc_hst_merged  = (hbt.incr2cum(n_hst_merged, DO_REVERSE=True)  * sarea / np.cos(subobslat_nh)).to('1').value
    n_cum_sc_aat         = (hbt.incr2cum(n_aat, DO_REVERSE=True) * sarea / np.cos(subobslat_nh)).to('1').value

# Now plot the line(s)

    if DO_AAT:
        plt.plot(r.to('micron').value, n_cum_sc_aat, label = 'q = ' + repr(q_i))

    if not(do_merged):
      if (DO_HST and DO_LAMBERT):
        plt.plot(r.to('micron').value, n_cum_sc_hst_lambert, label = 'q = ' + repr(q_i))

      if (DO_HST and DO_MIE):
        plt.plot(r.to('micron').value, n_cum_sc_hst_mie, label = 'q = ' + repr(q_i))

#    plt.plot(r.to('micron').value, n_cum_sc_aat, label = 'q = ' + repr(q_i))

    
# Plot the area between the lines

    if (DO_HST and DO_LAMBERT and DO_MIE):

# Plot a filled area between the Lambert and the Mie curves. This was a good idea, but in retrospect there are 
# more clear ways to do this.
# 
#      pl = polygon([r/um, reverse(r/um)], [n_cum_sc_hst_lambert, reverse(n_cum_sc_hst_mie)], /FILL_BACKGROUND, $
#                    FILL_COLOR=colors[i], /DATA, /CLIP, transpar=90); 8.2.1: /CLIP now works!

#      p = plot( (transpose([[r/um],[r/um]]))(*), (transpose([[n_cum_sc_hst_lambert], [n_cum_sc_hst_mie]]))[*], /OVER, $
#        color=colors[i], thick=10, transparency=90)			; 8.1: Do a fake 'polygon fill' by drawing a lot of lines

        if do_merged:
            plt.plot(r.to('micron'), n_cum_sc_hst_merged, label = 'q = {}'.format(q_i))

# Label r_c

plt.axvline(r_danger_1.to('micron').value,linestyle=':')

#plt.text(r_danger_1.to('micron').value * 1.3, ((p.yrange)[0])*10, '$r_c$', /DATA, /OVERPLOT)

# Draw dotted lines at the number of hits that SAS specified

plt.axhline(nhits_1, linestyle = ':')
plt.axhline(nhits_2, linestyle = ':')

# Done with plotting the individual lines... now finish and render the plot

plt.yscale('log')
plt.xscale('log')
plt.xlabel('Radius [$\mu$m]')
plt.ylabel('Total number N > r impacting NH during encounter')
plt.xlim((0.1, 1000))
plt.ylim(ylim)
plt.title(title)
plt.legend()
plt.show()

# Draw the word 'Danger' in the UR quadrant

#  te = text(r_danger_1 * 1.3/um, ((p.yrange)[1])/100, 'Danger', /DATA, /OVERPLOT)

# Draw the legend

#  pos_legend = [1,5d-3]        ; Upper-left corner of legend is put at this location
  
    
# Print HBT and date on image

#==============================================================================
# Now make a second plot, which is just of our various optical models for backscatter (HST).
# This is done only for one value of q (the smallest one in the array -- that is, last one calculated, first in array).
#==============================================================================

#  i++		; re-increment i (since loop was decrementing it!)
#  		; This puts i at the first index, which is the one for which this calculation is done.
#
# ; i = 0         ; For which value of q do we do? Set its index here.
#  
#  p = plot([3,4], xtitle = 'Radius [$\mu$m]', $
#      ytitle = 'Total number N > r impacting NH during encounter', $
#      xrange=[0.1,1000d],yrange=yrange, /XLOG, /YLOG, /NODATA, $
#      title = title)
#
#  po = polygon([r_danger_1/um,(p.xrange)[1],(p.xrange)[1],r_danger_1/um], $
#    [nhits_1, nhits_1,(p.yrange)[1],(p.yrange)[1]],/data, $
#    fill_color='pink',fill_transparency=0,/over, linestyle=' ')
#
#
#  p_arr = replicate(p, 3)
#  p_arr[0] = plot(r/um, n_cum_sc_hst_mie,     /OVER, name='q = ' + st(q[i]) + $
#    ', Mie, n = ' + st(n_refract) + ' + ' + st(abs(m_refract)) + 'i', color='blue',  thick=linethick[i], 
#    linestyle = '-')
#  p_arr[2] = plot(r/um, n_cum_sc_hst_lambert, /OVER, name='q = ' + st(q[i]) + $
#    ', Lambertian, albedo = ' + st(albedo), color='green', thick=linethick[i], linestyle = '__')
#  p_arr[1] = plot(r/um, n_cum_sc_hst_merged,  /OVER, name='q = ' + st(q[i]) + $
#    ', Mixed Lambertian + Mie, r$_{trans}$ = ' + st(r_trans/um) + ' $\mu$m', color='red',   
#    thick=linethick[i], linestyle = '-.')
#
#  pos_legend = [150,5d-3]        ; Upper-left corner of legend is put at this location
#                                 ; Except in 8.2.1, where it is the upper-right, ugh...
#                                 
#  l = legend(target=p_arr,position=pos_legend, /DATA)
#
## Draw the 'danger corner'
#
#  te = text(r_danger_1 / 1.3/um, ((p.yrange)[0])*10, '$r_c$', /DATA, /OVERPLOT)
#  pd = plot([1,1]*r_danger_1/um, [1d-20, 1d20],linestyle=':', /OVER)
#  te = text(r_danger_1 * 1.3/um, ((p.yrange)[1])/100, 'Danger', /DATA, /OVERPLOT)
#  pd = plot([1d-9, 1d9], [1,1]*nhits_1, linestyle = ':', /OVER)
#
#  fn = file_out + '_mie_vs_lambert_q' + st(q[i]) + '_bw.png' & p.save, fn & print, 'Wrote: ' + fn
#  fn = file_out + '_mie_vs_lambert_q' + st(q[i]) + '_bw.eps' & p.save, fn & print, 'Wrote: ' + fn
#  
#  stop
#
#
#end


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

# 
    
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