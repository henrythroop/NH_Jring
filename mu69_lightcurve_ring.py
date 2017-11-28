#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:36:00 2017

MU69_LIGHTCURVE_RING.PY -- Look at whether flat light curve could be caused by a ring.

I used these calculations in Keynote presentation 27-Nov-2017.
"MU69 ring lightcurve flatness HBT 27Nov17.key"

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

from   matplotlib.figure import Figurelsladfasdf

from   astropy.stats import sigma_clip

# HBT imports

import hbt

pi = math.pi

# OK, let's start by distributing a ring at 3000 km radius. 

a_ring_center = 1000*u.km     # Ring radius (center)

a_ring_center = 10000*u.km
da_ring = 500*u.km            # Ring width

# Define the ring. This seems like a really awkward way to do it wrt units. Is it not possible to do this better??

a_ring = np.array([(a_ring_center - da_ring/2).to(u.km).value, (a_ring_center + da_ring/2).to(u.km).value])*u.km

albedo_ring = 0.5
albedo_mu69 = 0.5

r_mu69 = 15*u.km

tau_ring = 3.5e-5

alam = 550*u.nm
d_hst = 2.4*u.m

dist_mu69 = 40*u.au

#ratio_flux_ring_mu69 = (2*pi*(a_ring[1]**2 - a_ring[0]**2) * tau_ring  /  (pi * r_mu69**2)).value

ratio_flux_ring_mu69 = 1

tau_ring = (pi * r_mu69**2) * ratio_flux_ring_mu69 / (2*pi*(a_ring[1]**2 - a_ring[0]**2))

print("Ratio Ring / MU69 = {}".format(ratio_flux_ring_mu69))
print("Ring radius = {}. Ring width = {}.".format(a_ring_center, da_ring))
print(" -> Tau_ring = {}".format(tau_ring))

# Q: What is PSF width of HST at this distance? Alex says it is 3700 km.

ang_psf = (1.22*alam / d_hst).decompose().value  # For angles, keep as a scalar, not astropy units

km_psf_mu69 = ang_psf * dist_mu69.to('km')   # Result: 1700 km

# Q: Would this be detected by MU69 stellar occultations, on the 16" telescopes?
# A: No. We had per-pixel SNR of like 5. Shadow speed was __ km/sec.

# Q: Would this be detected on the Gemini telescopes?

# Q: Could this be detected on the HST PSF? 
# A: Having this much flux this far out would probably be pretty non-stellar. 
# A: But hang on. I can put this arbitrarily close inward. e.g., at tau = 1e-4, and 1000 km, width 500 km.

# Q: What is the I/F limit from the Gemini observations?
