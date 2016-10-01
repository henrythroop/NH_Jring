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
import numpy as np
     
file_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

print 'Reading: ' + file_save
 
lun = open(file_save, 'rb')
t = pickle.load(lun)
lun.close()

# Repair various columns. For these text columns, the deal is that I had made them too short
# and strings were getting chopped. So, these commands replace them with larger columns -
# same type, but longer strings.

col = t['x_pos_star_image'].astype('S20000')
t.replace_column('x_pos_star_image', col)

col = t['y_pos_star_image'].astype('S20000')
t.replace_column('y_pos_star_image', col)

col = t['bg_argument'].astype('S30')  # Expand argument column from str1 to str30
t.replace_column('bg_argument', col)

col = t['Comment'].astype('S100')    # expand comment column from str13 to str100
t.replace_column('Comment', col)

col = t['x_pos_star_cat'].astype('S10000')
t.replace_column('x_pos_star_cat', col)

col = t['y_pos_star_cat'].astype('S10000')
t.replace_column('y_pos_star_cat', col)

stop


groups = astropy.table.unique(t, keys=(['Desc']))['Desc']
        
# Set initial location to first file in first group

index_group = 7
index_image  = 36

groupmask = (t['Desc'] == groups[index_group])

t_group = t[groupmask]
tg = t_group[index_image]  # Grab this, read-only, since we use it a lot.

# Now fix the tables. We want to go thru all of the x_pos_star_cat, and if there is no closing ), then 
# turn off the is_navigated flag, and delete the x_pos_star_cat entry, or make it longer.

#t = t_group # We use this a lot, so make it shorter
for i in range(np.size(t)):
    field = t['x_pos_star_cat'][i]
    if ('(' in field) and (')' not in field):
        print "mismatch found in record {}; setting to invalid.".format(i)
        t['is_navigated'][i] = False  # Invalidate the navigation for this record

####

# Now write the file back out. Don't overwrite existing file.

file_save = file_save.replace('.pkl', '.tmp.pkl')

lun = open(file_save, 'wb')
pickle.dump(t, lun)
lun.close()

print 'Wrote: ' + file_save 
