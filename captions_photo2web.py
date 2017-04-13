#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:47:27 2017

@author: throop
"""

# Process photo2web captions files 

import os
import glob
from subprocess import call  # Best way to call a function from python
from subprocess import check_output # Call, and return the output

# 
# 

file_captions = 'captions.txt'
dir_originals = 'originals'

# With no arguments, take all of the files in the 'originals' dir, and one-by-one process them.
# Except if originals doesn't exist, then use the current directory.

if (os.path.isfile(dir_originals)):
    dir = os.path.curdir + '/' + dir_originals
else:
    dir = os.path.curdir

dir = '/Users/throop/photos/Trips/Test'

file = os.path.join(dir, 'originals', file_captions)

images = [glob.glob("*.jpg"), 
          glob.glob("*.jpeg"),
          glob.glob("*.JPG"),
          glob.glob("*.JPEG")]

images = [item for sublist in images for item in sublist]

# Read the existing captions file
# It is 2 * N lines long -- with newlines separating each entry

lun = open(os.path.join(dir, file_captions), 'r') # 'r' = read
captions = lun.readlines()
lun.close()

# Read all the captions from JPEG files

captions_new = []

for image in images:
    caption = check_output(['exiftool', '-Description', '/Users/throop/photos/Trips/Test/originals/001_IMG_5456.jpg'])
    
# .extend(glob.glob("*.JPG")).extend(glob.glob("*.jpeg")).extend(glob.glob("*.JPEG"))
    
# With multiple arguments:
#   Read the captions.txt
#   Read the list of files
#   Process the appropriate one (or more)
#   Write out the file
   
# First see if we have a 