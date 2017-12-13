#!/usr/bin/env python
################################################################################
# production_rate.py
#
# Syntax:
#   python production_rate.py --param=value ... normal_if.tab directories...
#
# This program generates the calibration factors "production_0", equal to the
# rate of dust production (number of particles per square km per second) at the
# surface of a moon in the Pluto system, for particles in bin 0 (centered on 1
# mm radius).
#
# Inputs:
#   normal_if.tab       replace by the tabulation of ring normal I/F vs. radius.
#   directories...      one or more directory trees that will be searched for
#                       header files and associated array and array2 files.
#
# Options:
#   --bodies=[...]      a body ID or list of body IDS for the source moon(s) of
#                       the ring being calibrated:
#                           a SPICE ID for styx, nix, kerberos, hydra;
#                           a small positive integer for one of the new moons,
#                               using charon-like properties (body_type = 2).
#                           a small negative integer for one of the new moons,
#                               using kbo-like properties (body_type = 1).
#                       Default is [902,903,904,905].
#   --ort=string        prefix to place in front of the names of new moons.
#                       Default 'ort22'.
#   --pattern=string    an optional regular expression that must appear inside
#                       the header's file path to qualify as a match.
#   --drag=n            the drag index value (0, 1, 2 for 10x, 3 for 100x, ...)
#                       to select for calibration.
#   --plot=''           an optional directory in which to save a PDF plot of the
#                       ring's radial profile; blank to suppress. Default is
#                       blank. If the plot is generated, the profile of normal
#                       I/F is also shown as a blank line.
#   --profile=name      optional name of another file tabulating ring normal I/F
#                       vs. radius. This is included on the plot as a thin black
#                       line.
#   --verbose=n         1 to list the header files of grids as they are loaded;
#                       0 to supress this information. Default 1.
#   --smooth=km         amount by which to smooth the calibration profile, in
#                       km. Default 1000.
#   --save=''           an optional directory in which to save tabulations of
#                       I/F for all the rings; blank to suppress. Default is
#                       blank.
#   --prefix=string     optional prefix to prepend to all output files. Default
#                       is "bigfour".
#
# The factors are described by a dictionary of dictionaries, indexed as follows:
#   production_0[case][speed][rho][q]
# where
#   case = 1, 2, 3 or 4
#   speed = 1 or 2
#   rho = 10.**(-1./3.), 1., or 10.**(1./3.)
#   q = 2.5, 3. or 3.5
#
# The structure is saved in a pickle file "production_0_dragN.pickle", where
# "N" is replaced by the drag factor used
#
# Mark Showalter, mshowalter@seti.org
# Version 1.0: Dec 23, 2014
# Version 1.1: Feb 17, 2015: Revised to handle match strings and drag selection.
# Version 1.2: Apr 3, 2015: Finalized plot options.
################################################################################

import sys
import numpy as np
import pylab
import pickle

from tabulation import Tabulation
from nhdust import *

S_STEP = 10.**(1./3.)
LOG_S_STEP = np.log(S_STEP)

BODY_NAMES = ['styx', 'nix', 'kerberos', 'hydra']
BODY_IDS = [905, 902, 904, 903]
BODY_NAME_DICT = {905: 'styx', 902: 'nix', 904: 'kerberos', 903: 'hydra'}
BODY_TYPES = ['nominal', 'kbo', 'charon']   # 0, 1, 2

RHO_VALUES = [1./S_STEP, 1., S_STEP]  # g/cm**3
Q_VALUES = [2.5, 3.0, 3.5]

PHOTBIN_MIN = -6
PHOTBIN_MAX =  3

# Parameters of radial profiles
DR = 100.
SAMPLES = 1000

################################################################################
# Interpret command arguments
################################################################################

bodies = [905, 902, 904, 903]
ort = 'ort22'
pattern = None
drag = 1
plot = ''
save = ''
profile = None
verbose = 1
smooth = 1000
prefix = 'bigfour'

args = []

for arg in sys.argv[1:]:
    if arg.startswith('--bodies='):
        exec(arg[2:])
    elif arg.startswith('--ort='):
        ort = arg[len('--ort='):]
    elif arg.startswith('--pattern='):
        pattern = arg[len('--pattern='):]
    elif arg.startswith('--drag='):
        exec(arg[2:])
    elif arg.startswith('--plot='):
        plot = arg[len('--plot='):]
    elif arg.startswith('--save='):
        save = arg[len('--save='):]
    elif arg.startswith('--profile='):
        profile = arg[len('--profile='):]
    elif arg.startswith('--verbose='):
        exec(arg[2:])
    elif arg.startswith('--smooth='):
        exec(arg[2:])
    elif arg.startswith('--prefix='):
        prefix = arg[len('--prefix='):]
    elif arg.startswith('--'):
        raise ValueError('Unrecognized option: ' + arg)
    else:
        args.append(arg)

if len(args) < 2:
    raise IOError('missing arg(s)')

body_ids = []
body_names = []
body_types = []
type_name = 'nominal'

for b in bodies:
    body_id = abs(b)
    body_name = ort + '-' + ('%04d' % body_id)

    if b > 900:
        body_type = 0
        body_name = BODY_NAME_DICT[body_id]
    elif b > 0:
        body_type = 2
        type_name = BODY_TYPES[body_type]
    else:
        body_type = 1
        type_name = BODY_TYPES[body_type]

    body_ids.append(body_id)
    body_names.append(body_name)
    body_types.append(body_type)

verbose = (verbose != 0)

plot_dir = plot
plot = (plot != '')

save_dir = save
save = (save != '')

smooth = float(smooth)

################################################################################
# Load ring profile
################################################################################

normal_if_file = args[0]

a = np.fromfile(normal_if_file, sep=' ')
a = a.reshape(a.size/2,2)
profile1_x = a[:,0]
profile1_y = a[:,1]

normal_if = Tabulation(profile1_x, profile1_y)
# normal_if(r) returns the normal I/F at r (in km)

if profile:
    a = np.fromfile(profile, sep=' ')
    a = a.reshape(a.size/2,2)
    profile2_x = a[:,0]
    profile2_y = a[:,1]

################################################################################
# Load the grids
################################################################################

big_dict = load_tree(args[1:], verbose=verbose, drag=drag, match=pattern)

# indexing is [body_id][body_type][case][speed][beta]

################################################################################
# Determine IF_norm/production_0 for each case, speed, rho, and q
################################################################################

if_norm_over_production_0_1d = {}
if_norm_over_production_0_2d = {}

for case in (1,2,3,4):
  if_norm_over_production_0_1d[case] = {}
  if_norm_over_production_0_2d[case] = {}

  for speed in (1,2):
    if_norm_over_production_0_1d[case][speed] = {}
    if_norm_over_production_0_2d[case][speed] = {}

    for irho in range(len(RHO_VALUES)):
      rho = RHO_VALUES[irho]
      if_norm_over_production_0_1d[case][speed][rho] = {}
      if_norm_over_production_0_2d[case][speed][rho] = {}

      for q in Q_VALUES:
        print 'case, speed, irho, q = ', case, speed, irho, q

        images = []
        profiles = []
        for ibody in range(len(body_ids)):
          body_id = body_ids[ibody]
          body_type = body_types[ibody]

          beta_dicts = big_dict[body_id][body_type][case][speed]
            # beta_dict[beta] = header_dict

          # Calibrate the 2-D grids to normal I/F
          headers_with_grids = reindex_betas(beta_dicts, rho, PHOTBIN_MIN,
                                                              PHOTBIN_MAX)
          image = if_norm_2d(headers_with_grids, q)
          profile = profile_from_grid(image, headers_with_grids[0],
                                             DR, SAMPLES, smooth)

          images.append(image)
          profiles.append(profile)

        if_norm_over_production_0_2d[case][speed][rho][q] = np.sum(images,
                                                                   axis=0)
        if_norm_over_production_0_1d[case][speed][rho][q] = profiles

################################################################################
# Derive production_0 dictionary and save it
################################################################################

radii = np.arange(SAMPLES) * DR
mask = (radii >= normal_if.domain()[0]) & (radii <= normal_if.domain()[1])
reference = normal_if(radii[mask])

production_0 = {}
for case in (1,2,3,4):
  production_0[case] = {}

  for speed in (1,2):
    production_0[case][speed] = {}

    for rho in RHO_VALUES:
      production_0[case][speed][rho] = {}

      for q in Q_VALUES:
        profiles = if_norm_over_production_0_1d[case][speed][rho][q]

        summed = np.sum(profiles, axis=0)

        ratio = reference / summed[mask]
        production_0[case][speed][rho][q] = ratio.min()

filename = prefix + ('_production_0_%s_drag%d' % (type_name, drag))

f = open(filename + '.pickle', 'wb')
pickle.dump(production_0, f)
f.close()

################################################################################
# Make plots if necessary
################################################################################

limits = (radii >= 20.e3) & (radii <= 100.e3)

if plot or save:
  for case in (1,2,3,4):
    for speed in (1,2):
      for irho in range(len(RHO_VALUES)):
        rho = RHO_VALUES[irho]
        for q in Q_VALUES:
            profiles = if_norm_over_production_0_1d[case][speed][rho][q]
            scale = production_0[case][speed][rho][q]

            filename = (prefix + '_speed%d_q%3.1f_case%d_%s_rho%d_drag%d') % \
                              (speed, q, case, type_name, irho, drag)

            if plot:
              pylab.close()
              if profile is not None:
                  pylab.plot(profile2_x/1000., profile2_y, 'k', lw=0.5)

              pylab.plot(profile1_x/1000., profile1_y, 'k')

              summed = 0.
              for profile in profiles:
                  pylab.plot(radii[limits]/1000., profile[limits] * scale)
                  summed = summed + profile * scale

              pylab.plot(radii[limits]/1000., summed[limits], 'k', lw=2)

              pylab.xlim(20,100)
              pylab.xlabel('Radius (1000 km)')
              pylab.ylabel('Normal I/F')

              pylab.savefig(os.path.join(plot_dir,filename) + '.pdf')

            if save:
              summed = 0.
              scaled = []
              for profile in profiles:
                scaled.append(profile * scale)
                summed = summed + scaled[-1]

              f = open(os.path.join(save_dir,filename) + '.tab', 'w')
              for i in range(len(radii)):
                f.write('%5.0f ' % radii[i])
                f.write('%10.4e ' % summed[i])
                for j in range(len(scaled)):
                  f.write('%10.4e ' % scaled[j][i])
                f.write('\n')
              f.close()

  pylab.close()

################################################################################
