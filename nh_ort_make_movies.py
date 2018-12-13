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

# '/Users/throop/Data/MU69_Approach/throop/stacks/stack_KALR_MU69_Hazard_L4_2018344_MU69_Approach_n28_z6.png'

zoom = 4

files_png  = glob.glob(os.path.join(dir_in, f'*z{zoom}.png'))

files_fits = glob.glob(os.path.join(dir_in, f'*z{zoom}.fits'))

# Make a PNG movie

doy     = {}
img_haz = {}
wcs     = {}
reqids   = []

for file in files_fits:
    hdu = fits.open(file)
    img_i = hdu[2].data  # 1=image; 2=field; 3=differential
    # plt.imshow(img_i, origin='lower')
    wcs_i = WCS(file)
    # plot_img_wcs(img_i, wcs, width=50, do_show=False)
    hdu.close()
    # plt.show()
    
    file_short = file.split('/')[7]  # Longer than reqid -- it's the full string
    reqid_i = '_'.join(file_short.split('_')[3:6])
    doy_i = file_short.split('_')[5]
    reqids.append(reqid_i)
    img_haz[reqid_i] = img_i
    doy[reqid_i] = doy_i 
    wcs[reqid_i] = wcs_i

s = sorted(doy.items(), key=operator.itemgetter(1))

reqids_sorted = []
for s_i in s:
    reqids_sorted.append(s_i[0])
    
reqids = reqids_sorted

# for reqid_i in reqids:
#     plot_img_wcs(img_haz[reqid_i], wcs[reqid_i], width=50, do_show=False, title=reqid_i)
#     plt.show()



# Start the movie
    
#%%
f = plt.figure(frameon=False, figsize=(4, 5), dpi=100)
canvas_width, canvas_height = f.canvas.get_width_height()
ax = f.add_axes([0, 0, 1, 1])
ax.axis('off')
fps = 5
width=100
file_out_mov = f'movie_f{fps}_w{width}.mp4'

cmdstring = ('ffmpeg', 
    '-y', '-r', f'{fps}', # overwrite, fps
    '-s', '%dx%d' % (canvas_width, canvas_height), # size of image string
    '-pix_fmt', 'argb', # format
    '-f', 'rawvideo',  '-i', '-', # tell ffmpeg to expect raw video from the pipe
    '-vcodec', 'mpeg4', file_out_mov) # output encoding
p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)

for reqid_i in reqids:
    
    # draw the frame
    
    plot_img_wcs(img_haz[reqid_i], wcs[reqid_i], width=width, do_show=False, title=reqid_i)
    plt.draw()

    # extract the image as an ARGB string
    string = f.canvas.tostring_argb()

    # write to pipe
    p.stdin.write(string)

# Finish up
p.communicate()

    
    