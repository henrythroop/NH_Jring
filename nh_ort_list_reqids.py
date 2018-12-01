#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 13:44:05 2018

@author: throop
"""

import os
import glob

def nh_ort_list_reqids():
    
    """
    This is a q&d function to list all of the reqids of NH KEM images.
    It lists the ReqID, and then how many WCS'd and Backplaned files are on the local disk.
    
    This file can be run easily from the shell.
    """
    
    dir_data = '/Users/throop/Data/MU69_Approach'
    
    dir_wcs        = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
    dir_backplanes = os.path.join(dir_data, name_ort, 'porter')
    
    dirs_wcs = glob.glob(os.path.join(dir_wcs, '*'))
    dirs_wcs = sorted(dirs_wcs)
    
    for dir in dirs_wcs:
        reqid = dir.split('/')[-1]
        list_wcs = glob.glob(os.path.join(dir, '*.fits'))
        num_wcs = len(list_wcs)
        
        dir2=dir.replace(dir_wcs, dir_backplanes)
        list_backplanes = glob.glob(os.path.join(dir2, '*.fits'))
        num_backplanes = len(list_backplanes)
                                               
        print(f'{reqid:35}: {num_backplanes:3} / {num_wcs:3}')
        
if __name__ == '__main__':
    nh_ort_list_reqids()
    
    