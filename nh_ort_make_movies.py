#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 11:53:43 2018

@author: throop
"""

import glob
import math
import os.path
import os

import astropy
from   astropy.io import fits
import astropy.table
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library
import pickle # For load/save
import subprocess

from  astropy.wcs import WCS
import scipy
import copy

import operator

# HBT imports

import hbt

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular
from   get_radial_profile_backplane import get_radial_profile_backplane
from   plot_img_wcs import plot_img_wcs
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   scipy.optimize import curve_fit
from   wcs_translate_pix import wcs_translate_pix, wcs_zoom

# HBT imports

import hbt

dir_in = '/Users/throop/Data/MU69_Approach/throop/stacks'

hbt.unload_kernels_all()
sp.furnsh('kernels_kem_prime.tm')

# '/Users/throop/Data/MU69_Approach/throop/stacks/stack_KALR_MU69_Hazard_L4_2018344_MU69_Approach_n28_z6.png'

zoom = 4

files_png  = glob.glob(os.path.join(dir_in, f'*z{zoom}.png'))

files_fits = glob.glob(os.path.join(dir_in, f'*z{zoom}.fits'))

# Make a PNG movie

doy     = {}
img_haz = {}
img_field = {}
img_diff = {}
wcs     = {}
img_haz_med   = {}
reqids   = []

for file in files_fits:
    print(f'Reading FITS file: {file}')
    hdu = fits.open(file)
    # img_i = hdu[2].data  # 1=image; 2=field; 3=differential
    # plt.imshow(img_i, origin='lower')
    wcs_i = WCS(file)
    # plot_img_wcs(img_i, wcs, width=50, do_show=False)
    # plt.show()
    
    file_short = file.split('/')[7]  # Longer than reqid -- it's the full string
    reqid_i = '_'.join(file_short.split('_')[3:6])
    doy_i = file_short.split('_')[5]
    reqids.append(reqid_i)
    
    img_haz[reqid_i] = hdu[0].data
    img_field[reqid_i] = hdu[1].data
    img_diff[reqid_i] = hdu[2].data
    
    # Calc the median of this frame
    
    img_haz_med[reqid_i] = np.nanmedian(img_haz[reqid_i][img_haz[reqid_i] > 0] )

    if (10 < img_haz_med[reqid_i] < 20):
        img_haz[reqid_i] *= 29.967 / 19.967  # Balance the 20-sec exposures
        img_diff[reqid_i] *= 29.967 / 19.967

    if (5 < img_haz_med[reqid_i] < 10):
        img_haz[reqid_i] *= 29.967 / 9.967  # Balance the 10-sec exposures
        img_diff[reqid_i] *= 29.967 / 9.967  

    doy[reqid_i] = doy_i 
    wcs[reqid_i] = wcs_i

    hdu.close()

s = sorted(doy.items(), key=operator.itemgetter(1))

reqids_sorted = []
for s_i in s:
    reqids_sorted.append(s_i[0])
    
reqids = reqids_sorted


    
# Start the movie
    
#%%
fps = 8

width=150  # Width of the image, in LORRI pixels. Usually keep at 150. 200 starts to show hot pixels, which distract.

cmap = 'Greys_r'
hbt.fontsize(18)

# If requested, just use a subset of the frames (e.g., first 20)

num_files_max = None

num_files_max = 37  # Comment out or set to None to use all frames. There is a big jump to 38. 37 is good, though it 
                    # ends at K-14d.

if (num_files_max):
    reqids = reqids[0:num_files_max]

# Duplicate the final reqid a couple of times, just so the movie will pause

for i in range(5):
    reqids.append(reqids[-1])
    
dir_frames = '/Users/throop/Data/MU69_Approach/frames'

vmax = 100 # For vertical scaling
dpi  = 50  # Size of the output picture. 200 for cathy. 50 for web. 100 for high-res web.

file_out_base = os.path.join(dir_frames, f'movie_w{width}_d{dpi}_{cmap}')

file_out_gif = f'{file_out_base}_n{len(reqids)}_f{fps}.gif'

for i,reqid_i in enumerate(reqids):

    # Convert from DOY to calendar day

    file_out_frame = f'{file_out_base}_{i:03}.png'
    
    doy = reqid_i.split('_')[-1][-3:]
    
    et = sp.utc2et(f'2018::{doy} 12:00:00')
    utc = sp.timout(et, "Month DD, YYYY", 20)
    utc = utc.replace(' 0', ' ')
    
    # Draw the frame

    f = plt.figure(frameon=False, figsize=(10, 5), dpi=dpi)  # Change dpi to change output size
    # f.patch.set_facecolor('pink')
    canvas_width, canvas_height = f.canvas.get_width_height()
    ax = f.add_axes([0, 0, 1, 1])
    ax.axis('off')
    # ax.set_facecolor('pink')

    plt.subplot(1,2,1)
    ax = plt.gca()
    # ax.set_facecolor('pink')
    
    str_dt = f'K{int(doy)-366} days'
    plot_img_wcs(img_haz[reqid_i], wcs[reqid_i], width=width, do_show=False, title=str_dt,
                 do_stretch=False, vmin=5, vmax=vmax, cmap=cmap, do_inhibit_axes=True,
                 do_plot_fiducial=False,)

    # plt.subplot(1,3,2)    
    # plot_img_wcs(img_field[reqid_i], wcs[reqid_i], width=width, do_show=False, title=reqid_i,
    #              do_stretch=False, vmin=5, vmax=vmax, cmap=cmap, do_inhibit_axes=True)
                 

    plt.subplot(1,2,2)    
    plot_img_wcs(img_diff[reqid_i], wcs[reqid_i], width=width, do_show=False, title=utc,
                 do_stretch=False, vmin=-vmax*0.5, vmax=vmax, cmap=cmap, do_inhibit_axes=True,
                 do_plot_fiducial=False,)

    plt.tight_layout()
    # plt.gca().set_facecolor('yellow')
    
    plt.savefig(file_out_frame)
    plt.show()
    
    print(f'Wrote frame {i:03}: {file_out_frame}')

str_convert = f'convert -delay {int(100/fps)} {file_out_base}*.png {file_out_gif}'
os.system(str_convert)

# print(str_convert)

print(f'Wrote: {file_out_gif}')

str_rm = f'rm {file_out_base}*.png'
os.system(str_rm)
    