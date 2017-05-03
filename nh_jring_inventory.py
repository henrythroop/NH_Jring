# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:21:46 2016

@author: throop
"""

import hbt
import os
import pickle
import astropy
import spiceypy as sp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

##########
# NH_Jring_invengory
# Take an inventory of all of the NH Jupiter rings images.
# Make a PDF file showing all of the images, along with a bit of info.
# Also, make a text file with tabular data.
#
# HBT 14-Jun-2016
# HBT 2-May-2017 Updated to Python3, and for new directory structure.
#
##########

# Possible additions: Phase angle [DONE]
#                     Distance from Jupiter. 
#                     Time before/after C/A

filename_save = 'nh_jring_read_params_5711.pkl' # Filename to save parameters in
#filename_save = 'nh_jring_read_params_100.pkl' # Filename to save parameters in

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'

dir_out = '/Users/throop/Data/NH_Jring/out/'

file_out_pdf   = 'nh_jring_inventory.pdf'   # PDF output file which is created
file_out_txt   = file_out_pdf.replace('pdf', 'txt')

sp.furnsh(file_tm) # Commented out for testing while SPICE was recopying kernel pool.

# Check if there is a pickle save file found. If it is there, go ahead and read it.

if os.path.exists(filename_save):
    print("Loading file: " + filename_save)
#            self.load(verbose=False) # This will load self.t
    
    lun = open(filename_save, 'rb')
    t = pickle.load(lun)
#            self.t = t # All of the variables we need are in the 't' table.

    lun.close()
#            t = self.t

# Find and process all of the FITS header info for all of the image files
# 't' is our main data table, and has all the info for all the files.

else:

# Search for all files. Ideally I would search for all files that don't have 'opnav',
# but the way the filter is set up, it is easier to search for all that do. Result is the same.
    
    t = hbt.get_fits_info_from_files_lorri(dir_images, pattern='opnav')

num_files = np.size(t)

print('Read ' + repr(num_files) + ' files.')

# XXX Now shorten the list arbitrarily, just for testing purposes

#t = t[100:112] # XXX shorten it! 100 does not get to 4x4 territory
    
# Get a list of unique observation descriptions (i.e., 'groups')

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

num_groups = np.size(groups)
      
# Set initial location to first file in first group

index_group  = 0
index_image  = 0

num_rows = 3
num_cols = 4

# Set the starting row & column position

row = 1
col = 1

i = 1       # Current image #, starting at 1
i_on_page = 1  # Current image #, on this page
i_page    = 1  # Current page number, within this group
i_file = 0  # Index of this file, within the group
i_group = 0 # Index of this group, within the list of groups

lines_out = [''] # This array is the output text file. It will be written monolithically at the end. 
                 # Less chance of leaving open.

plt.rc('image', cmap='Greys_r')  # Use a non-inverted colormap: White = stars and ring, etc.


# Start the PDF

pp = PdfPages(dir_out + file_out_pdf)

fs = 5.5 # Font size for the PDF. 2.7 is very tiny. 7 is good but text is a bit too large per cell.

# Now loop over every image

for i_group,group in enumerate(groups):
    
    mask = t['Desc'] == group
    
    files = t[mask]
  
    dt  = 0
    t0  = 0

    header = " ----- " + group + " ----- "
    
    print
    print(header)
    lines_out.append('')
    lines_out.append(header)
    
    # Starting a new group, we start in UL corner on new page
    
    row     = 1
    col     = 1
    i_on_page  = 1
    i_page     = 1

#    print("Generating output...")
    
    for i_file,file in enumerate(files): # Loop over files in this particular group

      if (t0 ==0): # Calculate time since previous image
          dt_s = 0
      else:
          dt_s = float(file['ET']) - t0
      
      m, s = divmod(dt_s, 60)
      h, m = divmod(m, 60)
      dt_str = "{:3d}h {:2d}m {:2d}s".format(int(h), int(m), int(s))
      if (dt_s == 0): dt_str = '--'
      dt_str = dt_str.replace(' 0h', '').replace(' 0m', '')
      t0 = float(file['ET'])
      
      # Create a super-short version of the filename (cut out the ApID)
      
      file_trunc = file['Shortname'].replace('lor_', '').replace('_0x630_sci', '').\
        replace('_0x633_sci', '').replace('_opnav', '').replace('.fit', '')
      
      utc_trunc = sp.et2utc(sp.utc2et(file['UTC']),'C', 0)
      
      # Print a line of the table, to the screen and the file
      
      line =  "{:>3}/{:<3}: {:>3}, {:>3},  {},   {},   {},  {:6.3f},{:>12s},   {:.1f} deg, {:<9}".format(int(i), 
                                                     int(num_files), int(i_group), 
                                                     int(i_file), file_trunc, file['Format'], utc_trunc, 
                                                     file['Exptime'], (dt_str), file['Phase']*hbt.r2d, file['Target'])
      print(line)
      lines_out.append(line)
      
      arr = hbt.read_lorri(file['Filename'], bg_method = 'Polynomial', bg_argument = 4, frac_clip = 1)
      
      arr = hbt.remove_brightest(arr, 0.99)
      arr = -hbt.remove_brightest(-arr, 0.99)

      # Plot the image to the PDF

      p = plt.subplot(num_rows, num_cols,i_on_page) # Image 1       
      plt.imshow(arr, interpolation='none')
#      plt.tight_layout()              # Reduces whitespace, but in this case pushes captions off the edge of page!
      a = plt.gca()                    # Get the axes
      a.get_xaxis().set_visible(False) # Do not plot the axes
      a.get_yaxis().set_visible(False)
 
      n1 = int(file['N1'])  # NAXIS1: This will be either 1024 [1x1], or 255 [4x4]
      scalefac = n1  / 1024. # Generate a scaling factor. This is 1.0 for a normal image, 
                             # or 0.25 for a 4x4. Use it for placing text() properly.
        
      # Generate the text to put next to each image on the PDF
      
      label1 = "{}/{}: {}, {}, {} s; {}".format(int(i), int(num_files), file_trunc, file['Format'], 
                file['Exptime'], dt_str)
      label2 = r'{}, {:.1f}$^\circ$, Group {}, File {}'.format(utc_trunc,  file['Phase']*hbt.r2d, 
                 int(i_group), int(i_file)) 
      
      plt.text(0, -90*scalefac, label1, fontsize = fs)
      plt.text(0, -30*scalefac, label2, fontsize = fs)
      
      # Generate the text to be on the header of each page
      
      if (i_on_page == 1):
          plt.text(-250*scalefac, 200*scalefac, group + ', group ' + repr(i_group) + '         ' +  
                   " page " + repr(i_page), 
                   fontsize=fs*2, rotation=90)
           
#      print "Just plotted string: " + str
      
      # Now update the row/column for the next image
            
      i         += 1      
      i_on_page += 1
      
      col = np.mod(i_on_page-1, num_cols)+1   # Calc new column number. The +1 is since they go 1, 2, 3   not   0, 1, 2
      if (col == 1):
          row += 1
      if (row == num_rows+1): # If we have just filled up the current page, and are starting a new one

#          print "NEW PAGE! Page full."
                    
          fig = plt.gcf()
          pp.savefig(fig) # Close the page, and start a new one, with the same group
          fig.clf()

          i_on_page = 1   # Start the new images at UL corner of new page
          i_page += 1
          row    = 1
          col    = 1
#          print
    if (i_on_page != 1):      
        fig = plt.gcf()    # If we have finished the loop for the current group, start a new page. 
                           # (NB this will double-paginate for N=12 files.)
        pp.savefig(fig)
        fig.clf()
    
#    print "NEW PAGE! Done with group."
       
#      plt.imshow(arr)
#      plt.show()
      
#nh_jring_inventory()

# http://stackoverflow.com/questions/2252726/how-to-create-pdf-files-in-python

#pp.savefig(plt.gcf()) # This generates page 1

print("Writing pdf...")
pp.close()
print("Wrote: " + dir_out + file_out_pdf)

print("Writing txt...")
np.savetxt(dir_out + file_out_txt, np.array(lines_out), fmt='%s')
print("Wrote: " + dir_out + file_out_txt)

print('Processed ' + repr(np.size(t)) + ' files.')


# Now make some plots of various quantities vs. time
# For ring itself (excluding the Io 'monitoring' obs), the total is N = 517 obs.
# Including the erroneous 'monitoring' ones, N = 571.
# And then thre is also LORRI and anything at random. And Himalia ring of Cheng et al.

jca_utc = '2007 FEB 28 05:41:48' # Looked up from GV
jca_et  = sp.utc2et(jca_utc)

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

for i_group,group in enumerate(groups):

    plt.subplots(1, 3, figsize=(13,4))

    mask = t['Desc'] == group

    t_days = (t['ET'][mask] - jca_et) / 86400
    phase  = t['Phase'][mask] * hbt.r2d
    dist   = t['dxyz'][mask] / 71492 # Distance from Jup, km
    res = 0.3 / 1024 * hbt.d2r * t['dxyz'][mask] # Spatial resolution, km


    plt.subplot(1,3,1)
    plt.plot(t_days, phase, marker = 'o', linestyle='none', label=group)

    plt.xlabel('Days from Jupiter C/A')
    plt.ylabel('Phase Angle')
    plt.xlim((-5,4))
    plt.ylim((00,180))
    
    
    plt.subplot(1,3,2)
    plt.plot(t_days, dist, marker = 'o', linestyle='none')
    plt.title(group + ', n=' + repr(np.sum(mask)))
    plt.xlim((-5,4))
    plt.ylabel('Distance [RJ]')
    plt.xlabel('Days from Jupiter C/A')
    plt.ylim((20, 110))

    plt.subplot(1,3,3)

    plt.plot(t_days, res, marker = 'o', linestyle='none')
    plt.xlim((-5,4))
    plt.ylabel('Resolution [km]')
    plt.xlabel('Days from Jupiter C/A')
    plt.ylim((10, 110))
    plt.show()

    
plt.plot((t['ET'] - jca_et) / 86400, t['Exptime'],linestyle='none', marker = '+')
plt.xlabel('Days from Jupiter C/A')
plt.ylabel('Exptime [s]')
    
# arr = hbt.read_lorri(dir_images + '/lor_0035122328_0x633_sci_1.fit')
# arr = hbt.read_lorri(dir_images + '/lor_0035282790_0x633_sci_1.fit')

stop
