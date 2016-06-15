# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:21:46 2016

@author: throop
"""

import hbt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

##########
# NH_Jring_invengory
# Take an inventory of all of the NH Jupiter rings images.
# Make a PDF file showing all of the images, along with a bit of info.
# Also, make a text file with tabular data.
#
# HBT 14-Jun-2016
##########

# Possible additions: Phase angle. 
#                     Distance from Jupiter. 
#                     Time before/after C/A

#filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in
filename_save = 'nh_jring_read_params_100.pkl' # Filename to save parameters in

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

file_out_pdf   = 'nh_jring_inventory.pdf'   # PDF output file which is created
file_out_txt   = file_out_pdf.replace('pdf', 'txt')

# Check if there is a pickle save file found. If it is there, go ahead and read it.

if os.path.exists(filename_save):
    print "Loading file: " + filename_save
#            self.load(verbose=False) # This will load self.t
    
    lun = open(filename_save, 'rb')
    t = pickle.load(lun)
#            self.t = t # All of the variables we need are in the 't' table.

    lun.close()
#            t = self.t

# Find and process all of the FITS header info for all of the image files
# 't' is our main data table, and has all the info for all the files.

else:
    
    t = hbt.get_fits_info_from_files_lorri(dir_images)

print 'Read ' + repr(np.size(t)) + ' files.'

# Now shorten the list arbitrarily, just for testing purposes

#t = t[0:12] # XXX shorten it! 100 does not get to 4x4 territory
    
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
fs = 3 # Font size

i = 1       # Current image #, starting at 1
i_on_page = 1  # Current image #, on this page
i_page    = 1  # Current page number, within this group
i_file = 0  # Index of this file, within the group
i_group = 0 # Index of this group, within the list of groups

lines_out = [''] # This array is the output text file. It will be written monolithically at the end. Less chance of leaving open.

plt.rc('image', cmap='Greys_r')  # Use a non-inverted colormap: White = stars and ring, etc.

# Start the PDF

pp = PdfPages(file_out_pdf)

# Now loop over every image

for i_group,group in enumerate(groups):
    
    mask = t['Desc'] == group
    
    files = t[mask]
  
    dt  = 0
    t0  = 0

    header = " ----- " + group + " ----- "
    
    print
    print header
    lines_out.append('')
    lines_out.append(header)
    
    # Starting a new group, we start in UL corner on new page
    
    row     = 1
    col     = 1
    i_on_page  = 1
    i_page     = 1
    
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
      
      file_trunc = file['Shortname'].replace('lor_', '').replace('_0x630_sci', '').replace('_0x633_sci', '')
      
      # Print a line of the table, to the screen
      
      line =  "{:>3}, {:>3},  {},   {},   {},  {:6.3f},{:>12s},   {:<9}".format(int(i_group), int(i_file), file_trunc, file['Format'], file['UTC'], 
                                                     file['Exptime'], (dt_str), file['Target'])
      print line
      lines_out.append(line)
      
      arr = hbt.get_image_nh(file['Filename'], bg_method = 'Polynomial', bg_argument = 4, frac_clip = 1)
      
      arr = hbt.remove_brightest(arr, 0.99)
      arr = -hbt.remove_brightest(-arr, 0.99)

      # Plot the image to the PDF

      p = plt.subplot(num_rows, num_cols,i_on_page) # Image 1       
      plt.imshow(arr)
      plt.tight_layout()
      a = plt.gca()                    # Get the axes
      a.get_xaxis().set_visible(False) # Do not plot the axes
      a.get_yaxis().set_visible(False)
 
      n1 = int(file['N1'])  # NAXIS1: This will be either 1024 [1x1], or 255 [4x4]
      scalefac = n1  / 1024. # Generate a scaling factor. This is 1.0 for a normal image, or 0.25 for a 4x4. Use it for placing text() properly.
        
      # Generate the text to put next to each image on the PDF
      
      label1 = "{}, {}, {} s; {}".format(file_trunc, file['Format'], file['Exptime'], dt_str)
      label2 = "{}, Group {}, File {}".format(file['UTC'], int(i_group), int(i_file)) 
      plt.text(0, -80*scalefac, label1, fontsize = fs)
      plt.text(0, -30*scalefac, label2, fontsize = fs)
      
      if (i_on_page == 1):
          plt.text(-120*scalefac, 1000*scalefac, group + ", page " + repr(i_page), fontsize=fs*2, rotation=90)
           
#      print "Just plotted string: " + str
      
      # Now update the row/column for the next image
            
      i         += 1      
      i_on_page += 1
      
      col = np.mod(i_on_page-1, num_cols)+1   # Calc new column number. The +1 is because they go 1, 2, 3   not   0, 1, 2
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
        fig = plt.gcf()    # If we have finished the loop for the current group, start a new page. (NB this will double-paginate for N=12 files.)
        pp.savefig(fig)
        fig.clf()
    
#    print "NEW PAGE! Done with group."
       
#      plt.imshow(arr)
#      plt.show()
      
#nh_jring_inventory()

# http://stackoverflow.com/questions/2252726/how-to-create-pdf-files-in-python

#pp.savefig(plt.gcf()) # This generates page 1

pp.close()
print "Wrote: " + file_out_pdf

np.savetxt(file_out_txt, np.array(lines_out), fmt='%s')
print "Wrote: " + file_out_txt

print 'Processed ' + repr(np.size(t)) + ' files.'

# arr = hbt.get_image_nh(dir_images + '/lor_0035122328_0x633_sci_1.fit')
# arr = hbt.get_image_nh(dir_images + '/lor_0035282790_0x633_sci_1.fit')
