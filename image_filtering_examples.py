# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 09:54:52 2016

@author: throop
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 08:42:25 2014
 
Contains 2-dimensional filters for the spatial filtering of images. Some of
the code below is taken from psychopy [Peirce JW (2009) Generating stimuli
for neuroscience using PsychoPy. Front. Neuroinform. 2:10.
doi:10.3389/neuro.11.010.2008]. I had trouble installing this on my macs, so
I just pilfered the appropriate functions (sorry!).

http://www.srmathias.com/image-filtering/
 
@author: Sam
"""
 
import numpy as np
from skimage import exposure
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import misc
import hbt
import astropy
from   astropy.io import fits

 
def butter2d_lp(shape, f, n, pxd=1):
    """Designs an n-th order lowpass 2D Butterworth filter with cutoff
   frequency f. pxd defines the number of pixels per unit of frequency (e.g.,
   degrees of visual angle)."""
    pxd = float(pxd)
    rows, cols = shape
    x = np.linspace(-0.5, 0.5, cols)  * cols / pxd
    y = np.linspace(-0.5, 0.5, rows)  * rows / pxd
    radius = np.sqrt((x**2)[np.newaxis] + (y**2)[:, np.newaxis])
    filt = 1 / (1.0 + (radius / f)**(2*n))
    return filt
 
def butter2d_bp(shape, cutin, cutoff, n, pxd=1):
    """Designs an n-th order bandpass 2D Butterworth filter with cutin and
   cutoff frequencies. pxd defines the number of pixels per unit of frequency
   (e.g., degrees of visual angle)."""
    return butter2d_lp(shape,cutoff,n,pxd) - butter2d_lp(shape,cutin,n,pxd)
 
def butter2d_hp(shape, f, n, pxd=1):
    """Designs an n-th order highpass 2D Butterworth filter with cutin
   frequency f. pxd defines the number of pixels per unit of frequency (e.g.,
   degrees of visual angle)."""
    return 1. - butter2d_lp(shape, f, n, pxd)
 
def ideal2d_lp(shape, f, pxd=1):
    """Designs an ideal filter with cutoff frequency f. pxd defines the number
   of pixels per unit of frequency (e.g., degrees of visual angle)."""
    pxd = float(pxd)
    rows, cols = shape
    x = np.linspace(-0.5, 0.5, cols)  * cols / pxd
    y = np.linspace(-0.5, 0.5, rows)  * rows / pxd
    radius = np.sqrt((x**2)[np.newaxis] + (y**2)[:, np.newaxis])
    filt = np.ones(shape)
    filt[radius>f] = 0
    return filt
 
def ideal2d_bp(shape, cutin, cutoff, pxd=1):
    """Designs an ideal filter with cutin and cutoff frequencies. pxd defines
   the number of pixels per unit of frequency (e.g., degrees of visual
   angle)."""
    return ideal2d_lp(shape,cutoff,pxd) - ideal2d_lp(shape,cutin,pxd)
 
def ideal2d_hp(shape, f, n, pxd=1):
    """Designs an ideal filter with cutin frequency f. pxd defines the number
   of pixels per unit of frequency (e.g., degrees of visual angle)."""
    return 1. - ideal2d_lp(shape, f, n, pxd)
 
def bandpass(data, highpass, lowpass, n, pxd, eq='histogram'):
    """Designs then applies a 2D bandpass filter to the data array. If n is
   None, and ideal filter (with perfectly sharp transitions) is used
   instead."""
    fft = np.fft.fftshift(np.fft.fft2(data))
    if n:
        H = butter2d_bp(data.shape, highpass, lowpass, n, pxd)
    else:
        H = ideal2d_bp(data.shape, highpass, lowpass, pxd)
    fft_new = fft * H
    new_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_new)))    
    if eq == 'histogram':
        new_image = exposure.equalize_hist(new_image)
    return new_image
 
#def test():
#    """Test the filters."""
#    orig_image = misc.face()[:,:,0]
##    orig_image = mpimg.imread('NW066.jpg')[:,:,0] #comment to use stock image
##    orig_image = np.random.random_sample((1000,1000)) #use noise instead
#    fft_orig = np.fft.fftshift(np.fft.fft2(orig_image))
#    recon_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_orig)))
# 
#   
#    plt.subplot(431)
#    plt.title('Original image')
#    plt.imshow(orig_image)
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(432)
#    plt.title('FFT (log transformed)')
#    plt.imshow(np.log(np.abs(fft_orig)))
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(433)
#    plt.title('Reconstructed image')
#    plt.imshow(recon_image)
#    plt.gray()
#    plt.axis('off')
#   
#   #####
#   
#    filt = butter2d_lp(orig_image.shape, 0.2, 2, pxd=43)
#    fft_new = fft_orig * filt
#    new_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_new)))
#    new_image = exposure.equalize_hist(new_image)
#   
#    plt.subplot(434)
#    plt.title('Lowpass filter')
#    plt.imshow(filt)
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(435)
#    plt.title('FFT (log transformed)')
#    plt.imshow(np.log(np.abs(fft_new)))
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(436)
#    plt.title('Filtered image (histogram equalised)')
#    plt.imshow(new_image)
#    plt.gray()
#    plt.axis('off')
#
######
#   
#    filt = butter2d_hp(orig_image.shape, 0.2, 2, pxd=43)
#    fft_new = fft_orig * filt
#    new_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_new)))
#    new_image = exposure.equalize_hist(new_image)
#   
#    plt.subplot(437)
#    plt.title('Highpass filter')
#    plt.imshow(filt)
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(438)
#    plt.title('FFT')
#    plt.imshow(np.abs(fft_new))
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(439)
#    plt.title('Filtered image (histogram equalised)')
#    plt.imshow(new_image)
#    plt.gray()
#    plt.axis('off')
#   
#    filt = butter2d_bp(orig_image.shape, 1.50001, 1.50002, 2, pxd=43)
#    fft_new = fft_orig * filt
#    new_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_new)))
#    new_image = exposure.equalize_hist(new_image)
#   
#    plt.subplot(4,3,10)
#    plt.title('Bandpass filter')
#    plt.imshow(filt)
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(4,3,11)
#    plt.title('FFT')
#    plt.imshow(np.abs(fft_new))
#    plt.gray()
#    plt.axis('off')
#    plt.subplot(4,3,12)
#    plt.title('Filtered image (histogram equalised)')
#    plt.imshow(new_image)
#    plt.gray()
#    plt.axis('off')
#   
#    plt.show()
#   
#if __name__ == '__main__':
#    test()
    

#####
 
#orig_image = misc.face()[:,:,0]
dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
file1 = 'lor_0034602723_0x630_sci_1.fit'

image_orig = hbt.remove_brightest(hbt.get_image_nh(dir_images + file1), 0.97, symmetric=True)

fft_orig = np.fft.fftshift(np.fft.fft2(image_orig))
recon_image = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_orig)))


# f,n = (0.1, 1) : removes a bit too much ring!
#       (0.1, 3) : pretty good. Just a bit of horizontal banding on bottom
#       (0.1, 10) - same as (0.1,3) *** BEST. Actually, has a bit too much 'ringing' artifacts. 0.1 3 better.
#       (0.1, 100) - starting to introduce low-frq horizontal banding artifacts
#       (0.5, 3) - lets thru way too much hf. Also ringing.
#       (0.1, 1) - too much ring removed
#       (0.01, 1) - fine but sfit is better. Shape just isn't a good match.
#       (0.01, 100) - fine but sfit is better


# power = 8: pretty great fit, really. The band at bottom is removed better than with the fft.
#         2: corners not removed well.
#         5: very nice removal, no doubt

f = 0.1 # 0.01 means only remove the v. low frq's
n = 2

power = 5
str_fft = ', f = {}, n = {}'.format(f,n)
str_sfit = ', power = {} '.format(power)

filt = butter2d_lp(image_orig.shape, f, n, pxd=43)
fft_new = fft_orig * filt
image_fft = np.abs(np.fft.ifft2(np.fft.ifftshift(fft_new)))

image_sfit = hbt.sfit(image_orig, power)

#new_image = exposure.equalize_hist(new_image)

fig, ax = plt.subplots(2, 2, figsize=(20, 13 ))
   
#plt.subplot(221)
#plt.title('FFT, Lowpass filter' + str_params)
#plt.imshow(filt)
#plt.gray()
#plt.axis('off')
#
#plt.subplot(222)
#plt.title('FFT (log transformed)')
#plt.imshow(np.log(np.abs(fft_new)))
#plt.gray()
#plt.axis('off')

min = np.min(image_orig-image_sfit)
max = np.max(image_orig-image_sfit)
r = (min,max)

plt.subplot(2,3,1)
plt.imshow(image_orig)
plt.gray()
plt.title('Original')
plt.axis('off')

plt.subplot(232)
plt.title('sfit' + str_sfit)
plt.imshow(image_sfit)
plt.gray()
plt.axis('off')

plt.subplot(233)
plt.title('Image - sfit' + str_sfit)
plt.imshow(image_orig - image_sfit, vmin=min, vmax=max)
plt.gray()
plt.axis('off')

plt.subplot(2,3,4)
plt.imshow(image_orig)
plt.gray()
plt.title('Original')
plt.axis('off')

plt.subplot(235)
plt.title('FFT, Low-pass' + str_fft)
plt.imshow(image_fft)
plt.gray()

plt.subplot(236)
plt.title('Image - FFT, low-pass' + str_fft)
plt.imshow(image_orig - image_fft, vmin=min, vmax=max)
plt.gray()
plt.axis('off')   