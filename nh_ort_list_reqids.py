#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 13:44:05 2018

@author: throop
"""

import os
import glob
import time

from astropy.table import Table
import numpy as np

def nh_ort_list_reqids():
    
    """
    This is a q&d function to list all of the reqids of NH KEM images.
    It lists the ReqID, and then how many WCS'd and Backplaned files are on the local disk.
    
    This file can be run easily from the shell.
    """
    
    dir_data = '/Users/throop/Data/MU69_Approach'
    
    dir_backplanes_s = os.path.join(dir_data, 'throop', 'backplaned')
    dir_backplanes_t = os.path.join(dir_data, 'throop', 'backplaned_tunacan')
    
    dir_wcs = os.path.join(dir_data, 'porter')
    
    dirs_wcs = glob.glob(os.path.join(dir_wcs, '*'))
    dirs_wcs = sorted(dirs_wcs)

    # Set up the output arrays
    
    num_wcs_arr          = []  # Number of files in the WCS dir
    num_backplanes_s_arr = []  # Number of sunflower backplanes for this reqid
    num_backplanes_t_arr = []  # Number of tunacan backplanes for this reqid
    reqid_arr            = []  # REQID 
    t_latest_arr         = []  # Time of newest WCS'd file, in UTC or local time
    ctime_arr            = []  # Time of newest WCS'd file, in Unix float format
    
    # Loop over every directory, and get # of files, latest file, etc.
    
    for dir in dirs_wcs:

        reqid_arr.append(dir.split('/')[-1])
        
        list_wcs = glob.glob(os.path.join(dir, '*.fits'))
        
        if list_wcs:
            num_wcs_arr.append(len(list_wcs))
            
            latest_file = max(list_wcs, key=os.path.getctime)   # This is a bit of python magic I didn't know of!
    
            ctime = os.path.getctime(latest_file)  
            ctime_arr.append(ctime)
            t_latest_arr.append(time.ctime(ctime))
            
            dir2=dir.replace(dir_wcs, dir_backplanes_s)
            list_backplanes_s = glob.glob(os.path.join(dir2, '*.fits'))
    
            dir2=dir.replace(dir_wcs, dir_backplanes_t)
            list_backplanes_t = glob.glob(os.path.join(dir2, '*.fits'))
    
            num_backplanes_t_arr.append(len(list_backplanes_t))
            num_backplanes_s_arr.append(len(list_backplanes_s))
        else:
            num_backplanes_s_arr.append(0)
            num_backplanes_t_arr.append(0)
            t_latest_arr.append(0)
            ctime_arr.append(0)
            num_wcs_arr.append(0)
                                               
    t = Table([reqid_arr, num_wcs_arr, num_backplanes_s_arr, num_backplanes_t_arr, t_latest_arr, ctime_arr], 
              names=('ReqID', '# WCS', '# Sun', '# Tuna', 'WCS Time', 'ctime'), 
              )
    
    # Add column for number remaining to process
    
    t['# UnSun'] = t['# WCS'] - t['# Sun']
    t['# UnTuna'] = t['# WCS'] - t['# Tuna']
          
    # Sort the table, and prepare a copy for printing
          
    t2 = t.copy()
    t2.sort(['ctime'])
    t2.remove_column('ctime')
    
    # Print it
    
    print(t2)
  
    # Print it, again, all rows, just in case Python truncates it on us

    for row in t2:
        print(row.as_void())

    # Sum and print.
    # I wish I could then do t2.add_row() to add this to the main column, but there is a column mismatch
    
    print(t2.groups.aggregate(np.sum))
    
if __name__ == '__main__':

    nh_ort_list_reqids()
    
