#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 13:13:02 2018

@author: throop
"""

import os
import numpy as np

SPEED_KM_SEC = 14. # km/sec

# =============================================================================
# This is a test file from MRS. It basically flies thru our grids, but using just a straight path, rather than the 
# actual path. This is for ballpark checks.
#
# MRS / HBT ~5-Dec-2018
#
# MRS:
# "I had a thought about how to check the approximate I/Fs of a set of dust files _after_they are generated. The program is 
# attached. It sums up the density along each trajectory and converts it to a rough I/F. The answer is only approximate, 
# depending upon how much of the I/F comes from outside the set of size bins tabulated, and it also uses a hard-wired 
# spacecraft speed of 14 km/sec. The result of a call to all_dust_ifs() should be a dictionary in which all the values 
# range from (perhaps) slightly higher than 2e-7 down to several times lower. If there is a factor-of-eight error, it'll 
# show up. My run on Doug's old sun10k models, using the alternate trajectory, 
# did return a reasonable set of I/Fs, so I suspect the problem is at your end."
# =============================================================================

def dust_if(filename, albedo):
    """The approximate I/F along a spacecraft trajectory."""

    a = np.fromfile(filename, sep=' ')
    a = a.reshape(-1,14)
    radii_mm = a[0][1:]
    time_sec = a[1:,0]
    number_per_km3 = a[1:,1:]

    number_per_km2 = np.sum(number_per_km3, axis=0) * SPEED_KM_SEC
    fill_factor = np.pi * radii_mm**2 * 1.e-12 * number_per_km2
    i_over_f = albedo * np.sum(fill_factor)
    print("worked")

    return i_over_f

def all_dust_ifs(dir):
    """A dictionary of estimated I/Fs keyed by file name."""

    filenames = os.listdir(dir)
    if_dicts = {}
    for filename in filenames:
        if not filename.endswith('.dust'): continue

        albedo = float(filename.split('_pv')[1][:4])
        filepath = os.path.join(dir, filename)
        if_dicts[filename] = dust_if(filepath, albedo)

    return if_dicts


if (__name__ == '__main__'):
    
    dir = '/Users/throop/data/ORT5/throop/deliveries/dph-tunacan3.5kinc55/for_mehoke/'
    
d = all_dust_ifs(dir)

# d
# {'alternate_DPH-sun10k_y2.2_q2.5_pv0.05_rho0.10.dust': 6.869748801841756e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.05_rho0.22.dust': 1.4602166003418275e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.05_rho0.46.dust': 2.0939644365006483e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.05_rho1.00.dust': 2.3798633771851026e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.10_rho0.10.dust': 6.869755555804557e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.10_rho0.22.dust': 1.4602162677875603e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.10_rho0.46.dust': 2.0939602204340494e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.10_rho1.00.dust': 2.3798703162278971e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.30_rho0.10.dust': 4.0268917153761016e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.30_rho0.22.dust': 6.798240662350574e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.30_rho0.46.dust': 1.4822312254813247e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.30_rho1.00.dust': 2.0845573919608175e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.70_rho0.10.dust': 4.026900502654504e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.70_rho0.22.dust': 6.798236929098157e-08,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.70_rho0.46.dust': 1.4822315576422864e-07,
#  'alternate_DPH-sun10k_y2.2_q2.5_pv0.70_rho1.00.dust': 2.0845631768359713e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.05_rho0.10.dust': 1.673015086995119e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.05_rho0.22.dust': 1.8625878560641785e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.05_rho0.46.dust': 1.8981581626920738e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.05_rho1.00.dust': 1.9218612385598457e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.10_rho0.10.dust': 1.6730176633449693e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.10_rho0.22.dust': 1.8625906210025047e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.10_rho0.46.dust': 1.898154948495052e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.10_rho1.00.dust': 1.9218651169598445e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.30_rho0.10.dust': 1.5776866181855216e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.30_rho0.22.dust': 1.690616226052403e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.30_rho0.46.dust': 1.8349245830585066e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.30_rho1.00.dust': 1.9067172481869563e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.70_rho0.10.dust': 1.5776897303096485e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.70_rho0.22.dust': 1.6906166883080805e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.70_rho0.46.dust': 1.8349215009622244e-07,
#  'alternate_DPH-sun10k_y2.2_q3.5_pv0.70_rho1.00.dust': 1.9067211431903379e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.05_rho0.10.dust': 7.388785529745e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.05_rho0.22.dust': 1.6284276571387452e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.05_rho0.46.dust': 2.208696707133956e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.05_rho1.00.dust': 2.425934252877508e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.10_rho0.10.dust': 7.388793137897491e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.10_rho0.22.dust': 1.6284270789213183e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.10_rho0.46.dust': 2.2086920638226499e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.10_rho1.00.dust': 2.4259422881787467e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.30_rho0.10.dust': 4.208839462460662e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.30_rho0.22.dust': 7.311874685273831e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.30_rho0.46.dust': 1.6529781375706152e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.30_rho1.00.dust': 2.1987742355599134e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.70_rho0.10.dust': 4.208848735515194e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.70_rho0.22.dust': 7.31187129688424e-08,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.70_rho0.46.dust': 1.6529785832802845e-07,
#  'alternate_DPH-sun10k_y3.0_q2.5_pv0.70_rho1.00.dust': 2.1987802211545684e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.05_rho0.10.dust': 1.693822631298002e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.05_rho0.22.dust': 1.9095450763381058e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.05_rho0.46.dust': 1.9379415346932232e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.05_rho1.00.dust': 1.9580862118887832e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.10_rho0.10.dust': 1.6938251391004654e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.10_rho0.22.dust': 1.9095478229460146e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.10_rho0.46.dust': 1.9379382711061943e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.10_rho1.00.dust': 1.9580903853831708e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.30_rho0.10.dust': 1.590034961311437e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.30_rho0.22.dust': 1.7116426288354182e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.30_rho0.46.dust': 1.8811842812772022e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.30_rho1.00.dust': 1.9466800036937016e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.70_rho0.10.dust': 1.590038033423338e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.70_rho0.22.dust': 1.7116429832740773e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.70_rho0.46.dust': 1.8811812038762258e-07,
#  'alternate_DPH-sun10k_y3.0_q3.5_pv0.70_rho1.00.dust': 1.9466839930612316e-07}