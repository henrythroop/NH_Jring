# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:35:09 2016

@author: throop
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.misc import imread
import os
import numpy as np
import hbt

# Make an array of images on a PDF file.
# Based on http://stackoverflow.com/questions/2252726/how-to-create-pdf-files-in-python , near the end.

pp = PdfPages("out.pdf")

fs = 5
i = 1 # Current image #
num_rows = 2
num_cols = 2

im = hbt.dist_center(1001)

for row in hbt.frange(1,num_rows, num_rows):
    for col in hbt.frange(1, num_cols, num_cols):
        
        p = plt.subplot(num_rows, num_cols,i) # Image 1
       
        plt.imshow(im)
        a = plt.gca()
        a.get_xaxis().set_visible(False) # We don't need axis ticks
        a.get_yaxis().set_visible(False)
        
#        plotImage(files[1]) 
        str = '{:d}, X={:d}, Y={:d}'.format(int(i), int(col), int(row))
        plt.text(0, -20, str, fontsize = fs)
        plt.text(0, -60, str, fontsize = fs)

        if (i == 1):
            plt.text(00, -50, 'TEST CASE', fontsize=20)
            
        plt.tight_layout()
        print "Just plotted string: " + str
        i+=1

fig1 = plt.gcf()


pp.savefig(fig1) # This generates page 1

fig1.clf()

i=1

for row in hbt.frange(1,num_rows, num_rows):
    for col in hbt.frange(1, num_cols, num_cols):
        
        p = plt.subplot(num_rows, num_cols,i) # Image 1
       
        plt.imshow(im)
        a = plt.gca()
        a.get_xaxis().set_visible(False) # We don't need axis ticks
        a.get_yaxis().set_visible(False)
        
#        plotImage(files[1]) 
        str = '{:d}, X={:d}, Y={:d}'.format(int(i*100), int(col*100), int(row*100))
        plt.text(0, -20, str, fontsize = fs)
        plt.text(0, -60, str, fontsize = fs)

        if (i == 1):
            plt.text(00, -50, 'TEST CASE', fontsize=20)
            
        plt.tight_layout()
        print "Just plotted string: " + str
        i+=1

pp.savefig(plt.gcf()) # This generates page 2
        
pp.close()
