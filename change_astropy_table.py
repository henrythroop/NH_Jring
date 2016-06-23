# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 10:14:24 2016

@author: throop
"""


# Just a short code to change the datatype of an astropy table. I've done this a few times so now I have it here
# to remind me how.
#
# HBT 20-Jun-2016

import pickle
import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table  
     
file_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

print 'Reading: ' + file_save
 
lun = open(file_save, 'rb')
t = pickle.load(lun)
lun.close()

col = t['x_pos_star_image'].astype('S20000')
t.replace_column('x_pos_star_image', col)

col = t['y_pos_star_image'].astype('S20000')
t.replace_column('y_pos_star_image', col)

file_save = file_save.replace('.pkl', '.tmp.pkl')

lun = open(file_save, 'wb')
pickle.dump(t, lun)
lun.close()

print 'Wrote: ' + file_save 
