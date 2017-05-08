#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 13:26:01 2017

@author: throop
"""

def nh_jring_load_objectlist(file):
    
    t = Table.read(file, format = 'csv')

    return t
