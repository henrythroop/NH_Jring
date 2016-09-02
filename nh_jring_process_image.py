#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 12:34:07 2016

@author: throop
"""

##########
# Process Image -- do all necessary bg subtraction on current image
##########
                                                                                                                                                    
    def nh_jring_process_image(self):

        method = self.var_option_bg.get()

#        print "process_image()"
        
        frac, poly = 0, 0
        
        if (method == 'Previous'):
            file_prev = self.t_group['Filename'][self.index_image-1]
#            print "file =      " + filename
            print "file_prev = " + file_prev
            image_bg = hbt.read_lorri(file_prev, frac_clip = 1.0, bg_method = 'None', autozoom=True)
            image_fg = self.image_raw
            image = image_fg - image_bg

        if (method == 'Next'):
            file_next = self.t_group['Filename'][self.index_image+1]
            image_bg = hbt.read_lorri(file_next, frac_clip = 1.0, bg_method = 'None', autozoom=True)
            image_fg = self.image_raw
            image = image_fg - image_bg
            
        if (method == 'Median'): # XXX not working yet
            file_prev = self.t_group['Filename'][self.index_image-1]
            image_bg = hbt.read_lorri(file_prev, frac_clip = 1.0, bg_method = 'None')
            image_fg = self.image_raw
            image = image_fg - image_bg

        if (method == 'Polynomial'):
         power = self.entry_bg.get()
         image = self.image_raw - hbt.sfit(self.image_raw, power) # Look up the exponenent and apply it 
                                                
        if (method == 'Grp Num Frac Pow'):  # Specify to subtract a specified group#/image#, mult factor, and sfit power.
                                            # I thought this would be useful, but it turns out we usually need to subtract
                                            # a median of multiple images -- not just one -- so this is not very useful.
                                            # Plus, the best power is usually 5, and the best frac can be calc'd
                                            # with a linfit.
        
            vars = self.entry_bg.get().split(' ')
            
            if (np.size(vars) == 0): # If no args passed, just plot the image
                power = 0
                frac  = 0
                image = self.image_raw

            if (np.size(vars) == 1): # One variable: interpret as exponent
                power = float(vars[0])
                frac  = 0
                image = self.image_raw
                image_bg = hbt.sfit(image, power)
                image = image - image_bg
                
            if (np.size(vars) == 2): # Two variables: interpret as group num and file num
                (grp, num) = vars
                frac  = 1
                power = 0
                
            if (np.size(vars)) == 3: # Three variables: interpret as group num, file num, fraction
                (grp, num, frac) = vars
                power = 0
                
            if (np.size(vars) == 4): # Four variables: Group num, File num, Fraction, Exponent
                (grp, num, frac, power) = vars
                 
            if int(np.size(vars)) in [2,3,4]:
               
                grp = int(grp)
                num = int(num)
                frac = float(frac)
                power = int(power)
                
                print "group={}, num={}, frac={}".format(grp, num, frac)
#            print "Group = {}, num{}, Name = {}".format(name_group, num, name)

                name_group = self.groups[grp]
                groupmask = self.t['Desc'] == name_group
                group_tmp = self.t[groupmask]
                filename_bg = group_tmp['Filename'][num]
                                        
                image_fg = self.image_raw
                image_bg = hbt.read_lorri(filename_bg, frac_clip = 1, bg_method = 'None')
                
                image = image_fg - float(frac) * image_bg                
                image = image - hbt.sfit(image, power)
        
        if (method == 'String'):

# Parse a string like "6/112-6/129", or "129", or "6/114", or "124-129" or "6/123 - 129"
# As of 8-July-2016, this is the one I will generally use for most purposes.
# 'String' does this:
#   o Subtract the bg image made by combining the named frames
#   o Subtract a 5th order polynomial
#   o Filter out the extreme highs and lows
#   o Display it.    

            str = self.entry_bg.get()
            str2 = str.replace('-', ' ').replace('/', ' ')

            vars = str2.split(' ')

#            print "str = " + repr(str2)
#            print "vars = " + repr(vars)
            
            if (np.size(vars) == 0):
                image = self.image_raw
                self.image_processed = image
                return
                
            if (np.size(vars) == 1):
                stray = hbt.nh_get_straylight_median(self.index_group, [int(vars[0])])  # "122" -- assume current group
                
            if (np.size(vars) == 2):
                stray = hbt.nh_get_straylight_median(self.index_group, hbt.frange(int(vars[0]), int(vars[1])))  # "122-129" 
                                                                                            # -- assume current group
 
            if (np.size(vars) == 3):
                stray = hbt.nh_get_straylight_median(int(vars[0]), hbt.frange(vars[1], vars[2])) # "5/122 - 129"
                
            if (np.size(vars) == 4):
                stray = hbt.nh_get_straylight_median(int(vars[0]), hbt.frange(vars[1], vars[3])) # "6/122 - 6/129"

            image_proc = hbt.remove_sfit(self.image_raw, 5)
            (stray_norm,coeffs) = hbt.normalize_images(stray, image_proc) # Normalize the images to each other ()
            image_sub = image_proc - stray_norm
            
            image = image_sub
            
            image_bg = stray_norm
            
        if (method == 'None'):
            image = self.image_raw

        # Scale image by removing stars, CRs, etc. But be careful not to clamp too much: the whole point is that 
        # we are trying to make the ring bright, so don't reduce its peak!

        image =  hbt.remove_brightest( image, 0.97, symmetric=True)   # Remove positive stars
        
        # Save the image where we can ustse it later for extraction. 

        self.image_processed = image

# For analysis, and perhaps eventual integration into the GUI: plot the image,
# and the background that I remove. Plot to Python console, not the GUI.
        
        DO_DIAGNOSTIC = True
        
        if (DO_DIAGNOSTIC):
            plt.rcParams['figure.figsize'] = 12,6

            plt.subplot(1,3,1)
            plt.imshow(hbt.remove_brightest(hbt.remove_sfit(self.image_raw,5), 0.95, symmetric=True))
            plt.title('Original')
            
            plt.subplot(1,3,2)
            plt.imshow(hbt.remove_brightest(image_bg, 0.95, symmetric=True))
            plt.title('Background')
            
            plt.subplot(1,3,3)
            plt.imshow(image)
            plt.title('Result, med=' + np.)
            
            plt.show()
