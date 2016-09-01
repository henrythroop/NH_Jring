#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 15:53:57 2016

@author: throop
"""

### THS IS JUST SOME CODE I THOUGHT MIGHT BE USEFUL TO HAVE. IT IS HOW TO REPOINT WCS COORDS ###
# Originally I was using this in the J-ring analysis code. But it turned out that I just shifted the image
# in pixels instead.        

######## TEST WCS OFFSETTING ##########

# Now convert the dx / dy offset (from opnav) into an RA/Dec offset.
# Put the new correct center position into WCS, so plots from now on will be correct.

        do_offset_wcs = False
        
        if (do_offset_wcs):
            m = w.wcs.piximg_matrix
            w.wcs.crval[0] -= (dy_opnav * (m[0,0] + m[1,0])) # RA. I am sort of guessing about these params, but looks good.
            w.wcs.crval[1] -= (dx_opnav * (m[1,1] + m[0,1])) # Dec
    
            ax1 = figs.add_subplot(1,2,2, projection=w) # nrows, ncols, plotnum. Returns an 'axis'
            fig1 = plt.imshow(np.log(image_polyfit))
            # Get the x and y pixel positions from WCS. The final argument ('0') refers to row vs. column order
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)       
            
            plt.plot(points_phot[:,0], points_phot[:,1], marker='o', ls='None', fillstyle='none', color='red', 
                              markersize=12)
    
            x_stars,        y_stars        = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)      
            x_stars_abcorr, y_stars_abcorr = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
            x_ring_abcorr, y_ring_abcorr = w.wcs_world2pix(radec_ring_abcorr[:,0]*r2d, radec_ring_abcorr[:,1]*r2d, 0)
    
    
            plt.plot(x_stars, y_stars, marker='o', ls='None', color='lightgreen', mfc = 'None', mec = 'red', \
                              ms=12, mew=2, label='Cat')
            plt.plot(x_stars_abcorr, y_stars_abcorr, marker='o', ls='None', mfc = 'None', mec = 'green', \
                              ms=12, mew=2, label = 'Cat, abcorr')
            plt.plot(x_ring_abcorr, y_ring_abcorr, marker='o', ls='-', color='blue', mfc = 'None', mec = 'blue', 
                              ms=12, mew=2, 
                     label = 'Ring, LT+S')
                    
            plt.xlim((0,1000))
            plt.ylim((0,1000))
            plt.title('WCS, center = ' + repr(w.wcs.crval))