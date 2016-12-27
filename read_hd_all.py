#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 21:14:05 2016

@author: throop
"""

def read_hd_all(verbose=False):
    '''
    Read HD star catalog into an Astropy Table.
    '''
    
    file_hd = '/Users/throop/occult/hdoccul/catalogs/hd/catalog.dat'

#; ** HD catalog file format, from Vizier's ReadMe file. **
#;
#; --------------------------------------------------------------------------------
#;    Bytes Format  Units   Label    Explanations
#; --------------------------------------------------------------------------------
#;    1-  6  I6     ---     HD       [1/272150]+ Henry Draper Catalog (HD) number
#;    7- 18  A12    ---     DM       Durchmusterung identification (1)
#;   19- 20  I2     h       RAh      Hours RA, equinox B1900, epoch 1900.0
#;   21- 23  I3     0.1min  RAdm     Deci-Minutes RA, equinox B1900, epoch 1900.0
#;       24  A1     ---     DE-      Sign Dec, equinox B1900, epoch 1900.0
#;   25- 26  I2     deg     DEd      Degrees Dec, equinox B1900, epoch 1900.0
#;   27- 28  I2     arcmin  DEm      Minutes Dec, equinox B1900, epoch 1900.0
#;       29  I1     ---     q_Ptm    [0/1]? Code for Ptm: 0 = measured, 1 = value
#;                                         inferred from Ptg and spectral type
#;   30- 34  F5.2   mag     Ptm      ? Photovisual magnitude (2)
#;       35  A1     ---     n_Ptm    [C] 'C' if Ptm is combined value with Ptg
#;       36  I1     ---     q_Ptg    [0/1]? Code for Ptg: 0 = measured, 1 = value
#;                                         inferred from Ptm and spectral type
#;   37- 41  F5.2   mag     Ptg      ? Photographic magnitude (2)
#;       42  A1     ---     n_Ptg    [C] 'C' if Ptg is combined value for this
#;                                     entry and the following or preceding entry
#;   43- 45  A3     ---     SpType   Spectral type
#;   46- 47  A2     ---    Intensity [ 0-9B] Photographic intensity of spectrum (3)
#;       48  A1     ---     Remarks  [DEGMR*] Remarks, see note (4)

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy
    from astropy.table import Table
    import astropy.table as table
    from astropy.coordinates import SkyCoord
    import math
    import hbt
    
    #d2r = math.pi * 2 / 360
    
    # Set up empty lists for the output data
    
    num_hd   = [] # HD ID. 1 .. 272150
    mag_ptg  = [] # Photographic mag
    mag_ptm  = [] # Photometric mag
    ra       = [] # RA in radians (1950)
    dec      = [] # Dec in radians (1950)
    spt      = [] # Spectral type
    
    f = open(file_hd,'r')
    
    while True:
        l = f.readline()
        if (l == ''): break
        val_num_hd = int(l[0:6])
        val_ra_h = l[18:20]
        val_ra_dm = l[20:23]
        val_dec_sign = l[23]
        val_dec_deg  = (l[24:26])
        val_dec_min  = (l[26:28])
                         
        try:
            val_mag_ptm = float(l[29:34])
        except ValueError:
            val_mag_ptm = np.nan
        
        try:    
            val_mag_ptg = float(l[36:41])
        except ValueError:
            val_mag_ptg = np.nan
        
        val_spt = l[42:45]
        
        ra_str =  val_ra_h + 'h' + repr(0.10 * int(val_ra_dm)) + 'm'
        dec_str  =  val_dec_sign + val_dec_deg + 'd' + val_dec_min
        
        val_ra_rad = (int(val_ra_h) * 15 + (0.1 * int(val_ra_dm) * 15 / 60.) ) * hbt.d2r
        val_dec_rad = (int(val_dec_deg) + int(val_dec_min) / 60.) * hbt.d2r
        if ('-' in val_dec_sign):
            val_dec_rad += -1
    #        
    #    c = SkyCoord(ra_str, dec_str) # This is *really* slow! Total bottleneck.
    #    
    #    c.dec.degree
    #    c.ra.degree
    
        spt.append(val_spt)
    #    ra.append(c.ra.radian)
    #    dec.append(c.dec.radian)
        ra.append(val_ra_rad)
        dec.append(val_dec_rad)
        mag_ptm.append(val_mag_ptm)
        mag_ptg.append(val_mag_ptg)
        num_hd.append(val_num_hd)
        
        if ((np.mod(val_num_hd,100) == 0)) and (verbose):
            print("HD {}, RA = {}, dec = {}, spt = {}".format(val_num_hd, ra_str, dec_str, val_spt))
        
    #    num_hd.append(l[0:6])
    
    t = Table([num_hd, ra, dec, mag_ptg, mag_ptm, spt], names=['ID', 'RA', 'Dec', 'Ptg', 'Ptm', 'Type'])
    
    return(t)
