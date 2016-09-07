#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 13:24:29 2016

@author: throop
"""
# General python imports

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import pickle # For load/save
from mpl_toolkits.mplot3d import Axes3D
import astropy.visualization
from matplotlib import cm



# HBT imports

import hbt

#==============================================================================
# This is a test routine to call hbt.nh_process_image().
# This is just an easy wrapper, so I can debug that routine.
#==============================================================================

# This routine uses the same pickle save file as nh_jring_gui.
# The nh_jring_process_image() also uses this.

plt.set_cmap('Greys_r')

file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
    
lun = open(file_pickle, 'rb')
t = pickle.load(lun)
lun.close()

# Process the group names. Some of this is duplicated logic -- depends on how we want to use it.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

#    groupname = 'Jupiter ring - search for embedded moons'
#    groupnum = np.where(groupname == groups)[0][0]

index_image = 0  # Raw image, number
index_group = 7  # Raw image, group
method = 'String'   # Next, Previous, String, etc.
argument = '8-15'  # 33, or 8/23, or 10-20, etc.
    
groupmask = (t['Desc'] == groups[index_group])
t_group = t[groupmask]	

file = t_group[index_image]['Filename'] 

image = hbt.read_lorri(file)

image_proc = hbt.nh_jring_process_image(image, method, argument, index_group, index_image)

#stop

image_stray = hbt.nh_get_straylight_median(index_group, hbt.frange(8,15))

nx = np.shape(image)[0]
ny = np.shape(image)[1]

#image_stray = '/Users/throop/data/NH_Jring/out/straylight_median_g7_n8..15_sfit5,5.pkl'

xvals = range(nx)
yvals = range(ny)
(x,y) = np.meshgrid(xvals,yvals)

stretch = astropy.visualization.PercentileInterval(90)  # PI(90) scales array to 5th .. 95th %ile. 

image_s5 = hbt.remove_sfit(image,5)

plt.rcParams['figure.figsize'] = 5,5

#plt.subplot(2,2,1)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, stretch(image_stray), rstride=10, cstride=10, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()

#plt.subplot(2,2,2)
plt.imshow(stretch(image_stray))
plt.show()

#plt.subplot(2,2,3)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, stretch(image_s5), rstride=10, cstride=10, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()

#plt.subplot(2,2,4)
plt.imshow(stretch(image))
plt.show()
