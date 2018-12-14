#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 09:50:23 2018

@author: throop
"""

### Do some one-off testing to compares exposures and Field vs. OpaNav subtraction


file_haz = '/Users/throop/Data/MU69_Approach/throop/backplaned/' + \
             'KALR_MU69_OpNav_L4_2018345/lor_0406828338_0x633_pwcs2_backplaned.fits'
file_field = '/Users/throop/Data/MU69_Approach/throop/backplaned/' + \
             'K1LR_MU69ApprField_115d_L2_2017264/lor_0368314618_0x633_pwcs2_backplaned.fits'

lun = fits.open(file_haz)
exptime_haz = lun['PRIMARY'].header['EXPTIME']
arr_haz = lun['PRIMARY'].data
lun.close()

lun = fits.open(file_field)
exptime_field = lun['PRIMARY'].header['EXPTIME']
arr_field = lun['PRIMARY'].data
lun.close()

# hbt.figsize((10,8))
# h = img_haz[key]
# f = img_field[key]

h = arr_haz
f = np.roll(arr_field, (-42,0), axis=(0,1))

if (zoom == 4):
    hc = h[600:800, 600:800]
    fc = f[600:800, 600:800]

if (zoom == 1):
    hc = h[200:250, 200:250]
    fc = f[200:250, 200:250]

# exptime_haz   = np.median( stack_haz  [reqid_i].t['exptime'] )
# exptime_field = np.median( stack_field[reqid_i].t['exptime'] )

#%%%

# plt.hist(hc.ravel(), bins=50, range=(0,100), alpha=0.5, color='blue', label=f'Data, exptime={exptime_haz}')
# # plt.hist((fc * exptime_haz / exptime_field).ravel(), bins=50, range=(0,100), alpha=0.5, color='red', label='Field')
# plt.hist(fc.ravel(), bins=50, range=(0,100), alpha=0.5, color='red', label=f'Field, exptime={exptime_field}')
# plt.legend()
# plt.show()

#%%%

# hc = h - 7  # Haz frame (actually opnav, but close enough)
hc = h
fc = f   # Field frame

plt.set_cmap('plasma')
hbt.figsize((18,12))

vmin = 6
vmax = 30
plt.subplot(2,4,1)
plt.imshow(fc, vmin=vmin, vmax=vmax)
plt.title(f'Field, {exptime_field} sec')

plt.subplot(2,4,2)
plt.imshow(hc, vmin=vmin, vmax=vmax)
plt.title(f'Opnav, {exptime_haz} sec')

plt.subplot(2,4,3)
plt.imshow(hc-fc, vmin=vmin, vmax=vmax)
plt.title('Opnav-Field')

ratio = (exptime_haz/exptime_field) # This is messed up -- backwards??

plt.subplot(2,4,4)
plt.imshow(hc - fc*ratio, vmin=vmin, vmax=vmax)
plt.title(f'Opnav-{ratio:4.3}*Field')

#

plt.subplot(2,4,5)
plt.hist(fc.ravel(), bins=30, range=(-50,50))
plt.axvline(0)
plt.title(f'Field, {exptime_field} sec')

plt.subplot(2,4,6)
plt.hist(hc.ravel(), bins=30, range=(-50,50))
plt.axvline(0)
plt.title(f'Opnav, {exptime_haz} sec')

plt.subplot(2,4,7)
plt.hist((hc-fc).ravel(), bins=30, range=(-50,50))
plt.axvline(0)
plt.title('Opnav-Field')

plt.subplot(2,4,8)
plt.hist((hc-fc*ratio).ravel(), bins=30, range=(-50,50))
plt.axvline(0)
plt.title(f'Opnav-{ratio:4.3}*Field')

print(f'Field: {file_field}')
print(f'Opnav: {file_haz}')

plt.show()
