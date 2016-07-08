# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 21:51:55 2016

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

file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
dir_images =         '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
dir_backplanes =     '/Users/throop/data/NH_Jring/out/'

plt.rcParams['figure.figsize'] = 12, 12

lun = open(file_pickle, 'rb')
t = pickle.load(lun)
lun.close()

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

groupmask = (t['Desc'] == groups[index_group])
t_group = t[groupmask]

####

index_group = 6
index_image = hbt.frange(181,184).astype(int) # Which frame from the group do we extract?

image_stray = hbt.nh_get_straylight_median(6, hbt.frange(185,188)) # 

image_arr = np.zeros((1024, 1024, 4))
for i in range(4):
    image = hbt.get_image_nh(t_group[index_image[i]]['Filename'], autozoom=True)
    image_arr[:,:,i] = image

image_sum = np.sum(image_arr,axis=2)
image_sum_sfit = hbt.remove_sfit(image_sum, 5)

(image_stray_scale,junk) = hbt.normalize_images(image_stray, image_sum_sfit)

final = image_sum_sfit - image_stray_scale
final = hbt.remove_brightest(final, 0.95, symmetric=True)
final = hbt.remove_sfit(final,5)

ax=plt.subplot(1,3,2)
name = (t_group[index_image[0]]['Shortname'] + ' .. ' + 
        t_group[index_image[3]]['Shortname']).replace('lor_', '').replace('_0x633_sci_1.fit','')
plt.title(name)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.imshow(final)
#plt.show()
###

index_image = hbt.frange(177,180).astype(int) # Which frame from the group do we extract?
image_stray = hbt.nh_get_straylight_median(6, hbt.frange(185,188))

image_arr = np.zeros((1024, 1024, 4))
for i in range(4):
    image = hbt.get_image_nh(t_group[index_image[i]]['Filename'], autozoom=True)
    image_arr[:,:,i] = image

image_sum = np.sum(image_arr,axis=2)
image_sum_sfit = hbt.remove_sfit(image_sum, 5)

(image_stray_scale,junk) = hbt.normalize_images(image_stray, image_sum_sfit)

final = image_sum_sfit - image_stray_scale
final = hbt.remove_brightest(final, 0.95, symmetric=True)
final = hbt.remove_sfit(final,5)

ax = plt.subplot(1,3,1)
name = (t_group[index_image[0]]['Shortname'] + ' .. ' + 
        t_group[index_image[3]]['Shortname']).replace('lor_', '').replace('_0x633_sci_1.fit','')
plt.title(name)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.imshow(final)
#plt.show()

####

index_image = hbt.frange(185,188).astype(int) # Which frame from the group do we extract?
image_stray = hbt.nh_get_straylight_median(6, hbt.frange(173,176))

image_arr = np.zeros((1024, 1024, 4))
for i in range(4):
    image = hbt.get_image_nh(t_group[index_image[i]]['Filename'], autozoom=True)
    image_arr[:,:,i] = image

image_sum = np.sum(image_arr,axis=2)
image_sum_sfit = hbt.remove_sfit(image_sum, 5)

(image_stray_scale,junk) = hbt.normalize_images(image_stray, image_sum_sfit)

final = image_sum_sfit - image_stray_scale
final = hbt.remove_brightest(final, 0.95, symmetric=True)
final = hbt.remove_sfit(final,5)

ax = plt.subplot(1,3,3)
name = (t_group[index_image[0]]['Shortname'] + ' .. ' + 
        t_group[index_image[3]]['Shortname']).replace('lor_', '').replace('_0x633_sci_1.fit','')
plt.title(name)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.imshow(final)
plt.show()


im1 = get_image_nh(35108923)
