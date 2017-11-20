#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Mon Nov 20 09:38:43 2017

@author: throop
"""

import math       
import matplotlib.pyplot as plt # pyplot
import numpy as np
from   astropy import units as u           # Units library
from   pymiecoated import Mie

def scatter_lambert(ang_phase, albedo=1):
    
    """
    Return the phase function of a Lambertian scatterer. Phase function is properly normalized.
    
    Parameters
    -----
    
    ang_phase:
        Phase angle. Radians, or astropy units.
    albedo:
        Albedo. Optional.
        
    """
    
    pi = math.pi
    
    p11_lambert		 = 8/(3 * pi) * (np.sin(ang_phase) + (pi - np.sin(ang_phase)) * np.cos(ang_phase))
  
  			          # But hang on -- how is this phase function right since it doesn't include diffraction! Ugh...
  			          # The real phase function should really be half this, I think. This properly integrates to 
  			          # 2 over 0 .. 2 pi -- but what *should* integrate to 2 is the phase function including 
  			          # the diffraction spike, meaning these P11's in backscatter would go down.
  			          # OK, halving it. That will have the effect of making particles darker, harder to see, and 
  			          # increasing N for lambert.
  
    p11_lambert      *= 0.5
  			          
    p11_lambert      *= albedo
    
    return(p11_lambert)
    
def scatter_mie_ensemble(nm_refract, n, r, ang_phase, alam, do_plot=False):
    
    """ 
    Return the Mie scattering properties of an ensemble of particles.
    
    The size distribution may be any arbitrary n(r).
    
    The returned phase function is properly normalized s.t. \\int(0 .. pi) {P sin(alpha) d_alpha} = 2.
    
    Parameters
    ------

    nm_refract:
        Index of refraction. Complex.
    n:
        Particle number distribution.
    r: 
        Particle size distribution. Astropy units.
    ang_phase: 
        Scattering phase angle (array).
    alam: 
        Wavelength. Astropy units.
        
    Return values
    ----
    phase:
        Summed phase curve -- ie, P11 * n * pi * r^2 * qsca, summed. Array [num_angles]
    qsca:
        Scattering matrix. Array [num_radii]. Usually not needed.        
    p11_out:
        Phase curve. Not summed. Array [num_angles, num_radii]. Usually not needed.
        
    """
    
    num_r = len(n)
    num_ang = len(ang_phase)
    
    pi = math.pi
        
    k = 2*pi / alam
    
    qmie = np.zeros(num_r)
    qsca = np.zeros(num_r)
    qext = np.zeros(num_r)
    qbak = np.zeros(num_r)
    qabs = np.zeros(num_r)
    
    p11_mie  = np.zeros((num_r, num_ang))
    
    # Calc Q_ext *or* Q_sca, based on whether it's reflected or transmitted light
     
    x = (2*pi * r/alam).to('1').value
       
    print('Doing Mie code')
    
    # Mie code doesn't compute the phase function unless we ask it to, by passing dqv.
     
    for i,x_i in enumerate(x):
        mie = Mie(x=x_i, m=nm_refract)  # This is only for one x value, not an ensemble
        qext[i] = mie.qext()
        qsca[i] = mie.qsca()
        qbak[i] = mie.qb()  
    
        for j, ang_j in enumerate(ang_phase):
      
            (S1, S2)  = mie.S12(np.cos(pi*u.rad - ang_j)) # Looking at code, S12 returns tuple (S1, S2).
                                                             # For a sphere, S3 and S4 are zero.
                                                             # Argument to S12() is scattering angle theta, not phase
      
            sigma = pi * r[i]**2 * qsca[i]
      
       # Now convert from S1 and S2, to P11: Use p. 2 of http://nit.colorado.edu/atoc5560/week8.pdf
      
            p11_mie[i, j]  = 4 * pi / (k**2 * sigma) * ( (np.abs(S1))**2 + (np.abs(S2))**2) / 2
    
       
       
    if (do_plot):
        for i, x_i in enumerate(x):
            plt.plot(ang_phase.to('deg'), p11_mie[i, :], label = 'X = {:.1f}'.format(x_i))
        plt.title("N = {}".format(nm_refract))
        plt.yscale('log')
        plt.legend(loc='upper left')
        plt.xlabel('Phase Angle [deg]')
        plt.show()
    
    # Multiply by the size dist, and flatten the array and return just showing the angular dependence.

    terms = np.transpose(np.tile(pi * n * r**2 * qsca, (num_ang,1)))
    
    phase_out = np.sum(p11_mie * terms, axis=0)
    
    # Remove any units from this
    
    phase_out = phase_out.value

    # Now normalize it. We want \int (0 .. pi) {P(alpha) sin(alpha) d_alpha} = 2
    # The accuracy of the normalized result will depend on how many angular bins are used.
    
    d_ang = ang_phase - np.roll(ang_phase,1)
    d_ang[0] = 0
    tot = np.sum(phase_out * np.sin(ang_phase) * d_ang.value)
                                                
    phase_out = phase_out * 2 / tot
    
    # Return everything
    
    return(phase_out, p11_mie, qsca)


# =============================================================================
# Now do a test case
# =============================================================================

def test_scatter_mie_ensemble():
    
    import hbt
    import numpy as np
    
    from scatter_mie_ensemble import scatter_mie_ensemble
    
    pi = math.pi
    
    # Set up the scattering properties
    
    n_refract = 1.33
    m_refract = -0.001
    nm_refract = complex(n_refract,m_refract)
    
    # Define the wavelength 
    
    alam = 500*u.nm
    
    # Set up the size distribution
    
    r_min = 0.1
    r_max = 20
    num_r = 30
    
    q    = 3.5
    r    = hbt.frange(r_min, r_max, num_r, log=True)*u.micron
    n    = r.value**-q
    
    # Set up the angular distribution
    
    num_ang = 55
    
    ang_phase = hbt.frange(0, pi, num_ang)*u.rad # Phase angle.
    
# =============================================================================
#   Call the function
# =============================================================================
    
    (phase, p11_out, qsca) = scatter_mie_ensemble(nm_refract, n, r, ang_phase, alam, do_plot=True)
    
    # Now make a plot pair. One plot with a phase curve, and one with an n(r) distribution.
    
    plt.subplot(1,2,1)
    plt.plot(ang_phase.to('deg'), phase)
    plt.xlabel('Phase Angle [deg]')
    plt.title('q = {}'.format(q))
    plt.yscale('log')
    
    plt.subplot(1,2,2)
    plt.plot(r.to('micron'), n)
    plt.xlabel('Radius [micron]')
    plt.ylabel('Number')
    plt.yscale('log')
    plt.title('q={}'.format(q))
    plt.xscale('log')
    plt.show()
