# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:21:46 2016

@author: throop
"""

import hbt
import matplotlib.pyplot as plt

#def nh_jring_inventory():
    
filename_save = 'nh_jring_read_params.pkl' # Filename to save parameters in

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

# Start up SPICE

cspice.furnsh(file_tm)

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

<<<<<<< HEAD
print 'Read ' + repr(np.size(t)) + ' files.'

# Now shorten the list arbitrarily, just for testing purposes

#t = t[0:12] # XXX shorten it! 100 does not get to 4x4 territory
=======
t = t[0:12] # XXX shorten it!
>>>>>>> parent of 88fda31... J-ring Inventory processor now complete. Adding python code, and sample output. Works well.
    
# Get a list of unique observation descriptions (i.e., 'groups')

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

num_groups = np.size(groups)
num_files  = np.size(t)
      
# Set initial location to first file in first group

index_group  = 0
index_image  = 0

num_rows = 3
num_cols = 3

# Set the starting row & column position

row = 1
col = 1
fs = 3 # Font size

i = 1 # Current image #, starting at 1
i_page = 1  # Current image #, on this page

# Start the PDF

pp = PdfPages("out.pdf")

# Now loop over every image

for group in groups:
    
    mask = t['Desc'] == group
    
    files = t[mask]
  
    dt = 0
    t0  = 0
    print
    print " ----- " + group + " ----- "
    
    # Starting a new group, we start in UL corner on new page
    
    row     = 1
    col     = 1
    i_page  = 1
    for file in files: # Loop over files in this particular group

      if (t0 ==0): # Calculate time since previous image
          dt = 0
      else:
          dt = float(file['ET']) - t0
          
      t0 = float(file['ET'])
      
      file_trunc = file['Shortname'].replace('lor_', '').replace('_0x630_sci', '')
      
<<<<<<< HEAD
      file_trunc = file['Shortname'].replace('lor_', '').replace('_0x630_sci', '').replace('_0x633_sci', '')
      
      # Print a line of the table, to the screen
      
      line =  "{:>3}/{:<3}: {:>3}, {:>3},  {},   {},   {},  {:6.3f},{:>12s},   {:<9}".format(int(i), int(num_files), int(i_group), int(i_file), file_trunc, file['Format'], file['UTC'], 
                                                     file['Exptime'], (dt_str), file['Target'])
      print line
      lines_out.append(line)
=======
      # Print a line of the table
>>>>>>> parent of 88fda31... J-ring Inventory processor now complete. Adding python code, and sample output. Works well.
      
      print "{},   {},   {},   {:5.3f},   {:8.1f},   {:<9},   {}".format(file_trunc, file['Format'], file['UTC'], 
                                                     file['Exptime'], dt, file['Target'], file['Desc'])
    
    
      arr = hbt.get_image_nh(file['Filename'], bg_method = 'Polynomial', bg_argument = 4, frac_clip = 1)
      
      arr = hbt.remove_brightest(arr, 0.99)
      arr = -hbt.remove_brightest(-arr, 0.99)

      # Now send an image to the PDF

      p = plt.subplot(num_rows, num_cols,i_page) # Image 1
       
      plt.imshow(arr)
      a = plt.gca()
      a.get_xaxis().set_visible(False) # We don't need axis ticks
      a.get_yaxis().set_visible(False)
        
<<<<<<< HEAD
      # Generate the text to put next to each image on the PDF
      
      label1 = "{}/{:<3}:  {},  {},  {} s; {}".format(int(i), int(num_files), file_trunc, file['Format'], file['Exptime'], dt_str)
      label2 = "{},  Group {}, File {}".format(file['UTC'], int(i_group), int(i_file)) 
      plt.text(0, -80*scalefac, label1, fontsize = fs)
      plt.text(0, -30*scalefac, label2, fontsize = fs)
      
      if (i_on_page == 1):
          plt.text(-120*scalefac, 1000*scalefac, group + ", page " + repr(i_page), fontsize=fs*2, rotation=90)
           
#      print "Just plotted string: " + str
=======
      str = '{:d}, X={:d}, Y={:d}'.format(int(i), int(col), int(row))
      str2 = "{}, {}, {}".format(file['Exptime'], file['Format'], file['Target'])
      plt.text(0, -20, str, fontsize = fs)
      plt.text(0, -60, str2, fontsize = fs)

      if (i_page == 1):
          plt.text(00, -50, group, fontsize=10)
            
      plt.tight_layout()
      print "Just plotted string: " + str
>>>>>>> parent of 88fda31... J-ring Inventory processor now complete. Adding python code, and sample output. Works well.
      
      # Now update the row/column for the next image
            
      i      += 1      
      i_page += 1
      
      col = np.mod(i, num_cols)+1   # Calc new column number. The +1 is because they go 1, 2, 3   not   0, 1, 2
      if (col == 1):
          row += 1
      if (row == num_rows+1):
          pp.savefig(plt.gcf()) # This generates page 1
          i_page = 1
          row    = 1
          col    = 1
          print "NEW PAGE! Page full."
          print
    
    pp.savefig(plt.gcf())
    print "NEW PAGE! Done with group."
       
#      plt.imshow(arr)
#      plt.show()
      
#nh_jring_inventory()

# http://stackoverflow.com/questions/2252726/how-to-create-pdf-files-in-python

pp.savefig(plt.gcf()) # This generates page 1

pp.close()
<<<<<<< HEAD
print "Wrote: " + file_out_pdf

np.savetxt(file_out_txt, np.array(lines_out), fmt='%s')
print "Wrote: " + file_out_txt

print 'Processed ' + repr(np.size(t)) + ' files.'

# arr = hbt.get_image_nh(dir_images + '/lor_0035122328_0x633_sci_1.fit')
# arr = hbt.get_image_nh(dir_images + '/lor_0035282790_0x633_sci_1.fit')
=======
>>>>>>> parent of 88fda31... J-ring Inventory processor now complete. Adding python code, and sample output. Works well.
