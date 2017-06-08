"""
     Calculate the hazard to the NH spacecraft based on the I/F or tau
     This code generates the plots that go into my NH Ring Hazard whitepaper, Nov 2011.
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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import wcsaxes
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

do_mie 	   = True
do_lambert  = True
do_merged	   = True			# Mie and Lambert, mixed and cross-faded

do_hst = True

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

alpha_hst	= math.pi*u.rad - phase_hst	    # Alpha = scattering angle, in degrees
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

if (do_aat):
    tau_norm  = 0.0035    # This is the input optical depth. It is determined from the occultation, from AAT_hbt.pro (part 7).
    yrange = [1e-10,1e10]
    file_out = 'pluto_dust_limits_aat_v2_rmax' + repr(rmax.to('micron')) + '_bw'
    title = 'Pluto dust limits, 2006 AAT occultation, $\tau$ < 0.0035'

if (do_hst):
    iof_norm	= 3e-7
    tau_norm	= iof_norm
    yrange = np.array([1e-10,1e10]) / 1e4
    file_out = 'pluto_dust_limits_hst_v2_rmax' + repr(rmax.to('micron')) + '_bw'
    title = 'Pluto inner region, HST 2012 imaging, I/F < 3$\times 10^{-7}$, $r_{max}$ = ' + repr(rmax.to('micron')) + \
        ' $\mu$m'

subobslat_nh  = -49.07*u.deg                 # 2015 Mar 26, Pluto from NH

# Create the plot

q         = np.array([1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7])			# 9 elements long
colors    = np.array(['black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black'])
linestyle = np.array(['-', ':', '--', '-.', '-:', '__', '-', ':', '--'])
linethick = np.array([1, 1, 1, 1, 1, 1, 3, 3, 3])
  
#p = plot([3,4], xtitle = 'Radius [$\mu$m]' ,$
#      ytitle = 'Total number N > r impacting NH during encounter',xrange=[0.1,1000d],yrange=yrange, /XLOG, /YLOG, /NODATA, $
#      title = title)
   
#p_arr = replicate(p, sizex(q)) 

# Draw a dotted line at the two danger sizes that SAS specified at 3-Nov-2011 workshop

r    =  hbt.frange(rmin, rmax, num_r, log=True)*u.micron  # Astropy bug? When I run frange(), it drops the units.
qmie = np.zeros(num_r)
qsca = np.zeros(num_r)
qext = np.zeros(num_r)
qbak = np.zeros(num_r)
qabs = np.zeros(num_r)
P11  = np.zeros(num_r)


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
x = (2*math.pi * r/alam).value
   
print('Doing Mie code')

# Mie code doesn't compute the phase function unless we ask it to, by passing dqv.
 
if do_mie:
    for i,x_i in enumerate(x):  
      mie = Mie(x=x_i,m=nm_refract)  # This is only for one x value, not an ensemble
      qext[i] = mie.qext()
      qsca[i] = mie.qsca()
      qbak[i] = mie.qb()
      (S1, S2)  = mie.S12(np.cos(math.pi*u.rad - phase_hst)) # Looking at code, S12 returns tuple (S1, S2).
                                                             # For a sphere, S3 and S4 are zero.
      k = 1
      sigma = 1
      # Now convert from S1 and S2, to P11: Use p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf
      
      P11[i]  = 4 * math.pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2
      
#      mie_single, x, nm_refract, dqxt, dqsc, dqbk, dg, xs1, xs2, dph, dqv=dqv 
  							# UK Mie code. m_refract must be neg. x must be < 12,000.
                              # XXX I don't get this. Is X a scalar, or vector?
       
#    qext = dqxt		# Q_ext
#  qsca = dqsc		# Q_sca
#  qbak = dqbk		# Q_back (??)
#  p11_mie  = dph[*]	# Phase function. Computed only if angles dqv is passed.
 
p11_lambert		 = 8/(3 * math.pi) * (np.sin(phase_hst) + (math.pi - np.sin(phase_hst)) * np.cos(phase_hst)).value	
  
  			          # But hang on -- how is this phase function right since it doesn't include diffraction! Ugh...
  			          # The real phase function should really be half this, I think. This properly integrates to 
  			          # 2 over 0 .. 2 pi -- but what *should* integrate to 2 is the phase function including 
  			          # the diffraction spike, meaning these P11's in backscatter would go down.
  			          # OK, halving it. That will have the effect of making particles darker, harder to see, and 
  			          # increasing N for lambert.
  
p11_lambert       *= 0.5d
  			          
p11_lambert 	  *= albedo

qsca_mie		= qsca
qsca_lambert		= 2	# I think this is right. Scattering efficiency = 2 since it can scatter and diffract?

p11_merged		= (p11_mie  * frac_mie) + (p11_lambert  * frac_lambert)
qsca_merged		= (qsca_mie * frac_mie) + (qsca_lambert * frac_lambert)

print, 'Done UK Mie code'

  if (do_aat) then qmie = qext
  if (do_hst) then qmie = qsca

# Draw a polygon for the filled dangerous area

  po = polygon([r_danger_1/um,(p.xrange)[1],(p.xrange)[1],r_danger_1/um], $
    [nhits_1, nhits_1,(p.yrange)[1],(p.yrange)[1]],/data, $
    fill_color='pink',fill_transparency=0,/over, linestyle=' ')

  for i = sizex(q)-1, 0, -1  do begin		; Loop downward, so that the low q's are plotted on top for visibility
  
    n    = powerdist(r, 1d, q[i], 1d)
    
# Now normalize this to be number per cm2, in cross-section.
 
    n_hst_lambert = n * tau_norm / total(n * pi * r^2 * qsca_lambert * p11_lambert) 	; Added P11 20-Nov-2012
    n_hst_mie     = n * tau_norm / total(n * pi * r^2 * qsca_mie     * p11_mie) 	; Added P11 20-Nov-2012
    n_aat	  = n * tau_norm / total(n * pi * r^2 * qsca_mie) 

    n_hst_merged  = n * tau_norm / total(n * pi * r^2 * qsca_merged * p11_merged)

# Now convert into total hitting s/c

    n_cum_sc_hst_lambert = incr2cum(n_hst_lambert) * sarea / cos(subobslat_nh)
    n_cum_sc_hst_mie 	 = incr2cum(n_hst_mie)     * sarea / cos(subobslat_nh)
    n_cum_sc_hst_merged  = incr2cum(n_hst_merged)  * sarea / cos(subobslat_nh)
    n_cum_sc_aat 	 = incr2cum(n_aat)         * sarea / cos(subobslat_nh)

# Now plot the line(s)

    if (do_aat) then $
      p_arr[i] = plot(r/um, n_cum_sc_aat, /OVER, name='q = ' + st(q[i]), color=colors[i], thick=thick[i], $
        linestyl=linestyle[i])

    if (not(do_merged)) then begin
      if ((do_hst) and (do_lambert)) then $
        p_arr[i] = plot(r/um, n_cum_sc_hst_lambert, /OVER, name='q = ' + st(q[i]), color=colors[i], thick=thick[i], $
	linestyle=linestyle[i])

      if ((do_hst) and (do_mie)) then $
        p_arr[i] = plot(r/um, n_cum_sc_hst_mie, /OVER, name='q = ' + st(q[i]), color=colors[i], thick=thick[i], $
	linestyle=linestyle[i])
    end
# Plot the area between the lines

   if ((do_hst) and (do_lambert) and (do_mie)) then $

# Plot a filled area between the Lambert and the Mie curves. This was a good idea, but in retrospect there are 
# more clear ways to do this.
# 
#      pl = polygon([r/um, reverse(r/um)], [n_cum_sc_hst_lambert, reverse(n_cum_sc_hst_mie)], /FILL_BACKGROUND, $
#                    FILL_COLOR=colors[i], /DATA, /CLIP, transpar=90); 8.2.1: /CLIP now works!

#      p = plot( (transpose([[r/um],[r/um]]))(*), (transpose([[n_cum_sc_hst_lambert], [n_cum_sc_hst_mie]]))[*], /OVER, $
#        color=colors[i], thick=10, transparency=90)			; 8.1: Do a fake 'polygon fill' by drawing a lot of lines

    if (do_merged) then $
      p_arr[i] = plot(r/um, n_cum_sc_hst_merged, /OVER, name='q = ' + st(q[i]), color=colors[i], thick=linethick[i], $
        linestyle=linestyle[i])


# Write data out to a file

  if (do_write_txt_file) then begin
    if (do_hst and do_mie)    then write_dust_file, 'HstMie',     q[i], r, n_cum_sc_hst_mie
    if (do_hst and do_lambert)then write_dust_file, 'HstLambert', q[i], r, n_cum_sc_hst_lambert
    if (do_merged)            then write_dust_file, 'HstMixed'     , q[i], r, n_cum_sc_hst_merged
    
    if (do_aat)               then write_dust_file, 'AAT',         q[i], r, n_cum_sc_aat
  end

  end

# Label r_c

  te = text(r_danger_1 / 1.3/um, ((p.yrange)[0])*10, '$r_c$', /DATA, /OVERPLOT)
  
  pd = plot([1,1]*r_danger_1/um, [1d-20, 1d20],linestyle=':', /OVER)

# Draw the word 'Danger' in the UR quadrant

  te = text(r_danger_1 * 1.3/um, ((p.yrange)[1])/100, 'Danger', /DATA, /OVERPLOT)

# Draw dotted lines at the number of hits that SAS specified

  pd = plot([1d-9, 1d9], [1,1]*nhits_1, linestyle = ':', /OVER)
  
# Draw the legend

  pos_legend = [1,5d-3]        ; Upper-left corner of legend is put at this location
  
  l = legend(target=p_arr,position=pos_legend, /DATA)

# Print HBT and date on image

  do_plot_date = 0

  if (do_plot_date) then $
    t = text(0.78, 0.015, 'HBT 4-Nov-2011', /NORM)
   
  fn = file_out + '.png' & p.save, fn & print, 'Wrote: ' + fn
  fn = file_out + '.eps' & p.save, fn & print, 'Wrote: ' + fn


#==============================================================================
# Now make a second plot, which is just of our various optical models for backscatter (HST).
# This is done only for one value of q (the smallest one in the array -- that is, the last one calculated, first in array).
#==============================================================================

  i++		; re-increment i (since loop was decrementing it!)
  		; This puts i at the first index, which is the one for which this calculation is done.

 ; i = 0         ; For which value of q do we do? Set its index here.
  
  p = plot([3,4], xtitle = 'Radius [$\mu$m]', $
      ytitle = 'Total number N > r impacting NH during encounter', $
      xrange=[0.1,1000d],yrange=yrange, /XLOG, /YLOG, /NODATA, $
      title = title)

  po = polygon([r_danger_1/um,(p.xrange)[1],(p.xrange)[1],r_danger_1/um], $
    [nhits_1, nhits_1,(p.yrange)[1],(p.yrange)[1]],/data, $
    fill_color='pink',fill_transparency=0,/over, linestyle=' ')


  p_arr = replicate(p, 3)
  p_arr[0] = plot(r/um, n_cum_sc_hst_mie,     /OVER, name='q = ' + st(q[i]) + $
    ', Mie, n = ' + st(n_refract) + ' + ' + st(abs(m_refract)) + 'i', color='blue',  thick=linethick[i], linestyle = '-')
  p_arr[2] = plot(r/um, n_cum_sc_hst_lambert, /OVER, name='q = ' + st(q[i]) + $
    ', Lambertian, albedo = ' + st(albedo), color='green', thick=linethick[i], linestyle = '__')
  p_arr[1] = plot(r/um, n_cum_sc_hst_merged,  /OVER, name='q = ' + st(q[i]) + $
    ', Mixed Lambertian + Mie, r$_{trans}$ = ' + st(r_trans/um) + ' $\mu$m', color='red',   thick=linethick[i], linestyle = '-.')

  pos_legend = [150,5d-3]        ; Upper-left corner of legend is put at this location
                                 ; Except in 8.2.1, where it is the upper-right, ugh...
                                 
  l = legend(target=p_arr,position=pos_legend, /DATA)

# Draw the 'danger corner'

  te = text(r_danger_1 / 1.3/um, ((p.yrange)[0])*10, '$r_c$', /DATA, /OVERPLOT)
  pd = plot([1,1]*r_danger_1/um, [1d-20, 1d20],linestyle=':', /OVER)
  te = text(r_danger_1 * 1.3/um, ((p.yrange)[1])/100, 'Danger', /DATA, /OVERPLOT)
  pd = plot([1d-9, 1d9], [1,1]*nhits_1, linestyle = ':', /OVER)

  fn = file_out + '_mie_vs_lambert_q' + st(q[i]) + '_bw.png' & p.save, fn & print, 'Wrote: ' + fn
  fn = file_out + '_mie_vs_lambert_q' + st(q[i]) + '_bw.eps' & p.save, fn & print, 'Wrote: ' + fn
  
  stop


end


#==============================================================================
#  Verify normalization of phase functions P11, etc.
#==============================================================================

def other():

  theta_0 = 0
  theta_1 = pi
  num_theta = 100

  theta = frange(theta_0, theta_1, num_theta)
  dtheta = theta[1] - theta[0]

  p11 = cos(theta)  		; lambertian phase function = cos(theta). This is the *surface* phase function --
  				; not the phase function for the entire sphere, which needs to be convolved with the
				; crescent shape!

  p11 = 1			; Integral = 2 for this value (isotropic scattering)
  p11 = 8/(3d * pi) * (sin(theta) + (pi - theta) * cos(theta))			
  				; Integral = 2 for this expression (disk-integrated scattering from lambert sphere)
  				; Madhusudhan & Burrows, 2012, http://arxiv.org/pdf/1112.4476v1.pdf, eq33

  int = total(p11 * sin(theta) * dtheta)		; This integral should be 2, as per TE98 eq 16.
 
  stop

#==============================================================================
# Calculate the largest particle size in the ring.
#==============================================================================

# The philosophy here is that according to Gruen et al 1985, the mass flux of interplanetary dust takes a steep dropoff
# above 1 gram (cf fig. 3, GZF85, p. 10). So, I assume that a 1 gram particle is the largest that will hit the Pluto ring
# system. Maybe this is a bit ad hoc... their figure three shows a knee, but that may be just due to how it's plotted, and 
# the real power dist may continue to go on forever.
# Flux of 1 gram particles is 2.2d-15/m2/sec at 1 AU. Maybe it's down by 100x at Pluto (my guess).
# That's basically one per day hitting Charon.

  v_imp = 5 * km
  m_imp	= 1d * gram
  k_ej 	= 1d-8
  m_ej 	= 0.5 * m_imp * v_imp^2 * k_ej
  r0	= 0.01 * um
  r1 	= 100*cm
  q_ej  = 3.5
  num_r	= 100
  r	= frange(r0, r1, num_r, /log)
  n	= powerdist(r, m_ej, q_ej, rho)	; Create the ejecta. Yes, most of this will not stay in the ring, but we are
  						; conservative and ignore that.
  rho	= 1
  m	= total(4/3d * pi * r^3 * rho)

  n	= n * (m_ej / m)			; Normalize the ejecta number to the proper mass
  r_max	= r[max(where(n gt 1))]

  print, 'For impactor size ' + st(m_imp) + ' grams, the max ejecta size is ' + st(r_max/um) + ' um.'

end

;;;

def write_dust_file(name, q, r, n)
    
  
  file_out = 'DustDist' + name + '_q' + st(q) + '_rmax' + st(max(r)/um) + '.txt'

  openw, lun, file_out, /get_lun
  printf, lun, 'q = ' + st(q) + ', ' + name
  printf, lun, 'Radius_microns Cumulative_hits_N>r_during_encounter'

  for i = 0, sizex(n)-1 do printf, lun, r[i]/um, n[i]
  close, lun & free_lun, lun
  print, 'Wrote: ' + file_out
end
