#!/usr/bin/env python
################################################################################
# calibrate_3d.py
#
# Syntax
#   python calibrate_3d.py --param=value ... body_id destination directories...
#
# This program generates a set of 3D arrays, calibrated in units of particles
# per cubic km, for the six designated size bins. The file is saved as a NumPy
# save file '.npz' containing two arrays:
#   info = np.array([imax, jmax, kmax, xcell, ycell, zcell])
#   density = a 4-D array of particle density (number per km**3), indexed by
#       [bin - bin_min, i, j, k]
#
# Each file's name is name_speed*_q*_case*_[nominal|charon|kbo]_rho*.dust
#
# Inputs:
#   body_id             the body ID of the moon being calibrated:
#                           a SPICE ID for styx, nix, kerberos, hydra;
#                           a small integer for one of the new moons.
#   destination         a path to the destination directory, where the dust
#                       files will be saved.
#   directories...      one or more directories to be searched for header and
#                       array files.
#
# Options:
#   --prefix=string     prefix to place in front of the names of new moons.
#                       Default 'ort22'.
#   --pattern=string    an optional regular expression that must appear inside
#                       the header's file path to qualify as a match.
#   --drag=n            the drag index value (0, 1, 2 for 10x, 3 for 100x, ...)
#                       to select for calibration.
#   --production=name   name of production_0 file, default production_0.pickle.
#   --validate=flag     1 to validate the calibrations; 0 otherwise.
#   --profile=name      optional name of a ring profile by which to validate
#                       the calibrations and/or to appear on the plots.
#   --inner=flag        1 for an inner moon; 0 for an outer moon. Default 0.
#                       If inner = 1, then only case 1 input files are expected,
#                       but cases 2-4 will be used to scale the model.
#   --plot=dir          optional directory in which to save a PDF plot of each
#                       ring profile. The name will match that of the 'npz' file
#                       except for suffix '_profile.pdf'. If validate is
#                       not blannk, the validation profile will be added in red.
#                       If blank, then no plots are generated
#   --images=dir        an optional directory in which to save an image of each
#                       ring. If blank, then no images are generated.
#   --axis=axes         an integer, list ortuple containing one or more axis
#                       values:
#                           if 2 or -1 appears, a z-axis view;
#                           if 0, an x-axis view;
#                           if 1, a y-axis view.
#   --compare=name      an optional name of a ring's I/F tabulation. Used only
#                       if plot is nonzero. Adds a second reference plot in
#                       green.
#   --save=dir          optional directory in which to save a tabulation of each
#                       ring profile. The name will match that of the 'npz' file
#                       except for suffix '_profile.tab'. The file has two ASCII
#                       columns, radius in km and ring normal I/F. If blank,
#                       then no tables are generated
#
# Mark Showalter, mshowalter@seti.org
# Version 1.0: December 23, 2014
# Version 1.1 1/3/14: Revised to add destination arg and to use npz files, added
#   validation option.
# Version 1.2: Feb 17, 2015: Revised to handle match strings and drag selection;
#   added inner flag so that inner moons are calibrated according to all four
#   cases.
# Version 1.3: 4/3/15: Added plot and image options.
################################################################################

import sys
import os
import numpy as np
import pickle

from nhdust import *
from scipy.ndimage.filters import gaussian_filter
from tabulation import Tabulation

################################################################################
# Interpret command arguments
################################################################################

prefix = 'ort22'
pattern = None
drag = 1
production = 'production_0.pickle'
validate = 0
profile = ''
inner = 0
plot = ''
images = ''
axis = -1
compare = ''
save = ''
args = []

for arg in sys.argv[1:]:
    if arg.startswith('--prefix='):
        prefix = arg[len('--prefix='):]
    elif arg.startswith('--pattern='):
        pattern = arg[len('--pattern='):]
    elif arg.startswith('--production='):
        production = arg[len('--production='):]
    elif arg.startswith('--validate='):
        exec(arg[2:])
    elif arg.startswith('--profile='):
        profile = arg[len('--profile='):]
    elif arg.startswith('--drag='):
        exec(arg[2:])
    elif arg.startswith('--inner='):
        exec(arg[2:])
    elif arg.startswith('--plot='):
        plot = arg[len('--plot='):]
    elif arg.startswith('--images='):
        images = arg[len('--images='):]
    elif arg.startswith('--axis='):
        exec(arg[2:])
    elif arg.startswith('--compare='):
        compare = arg[len('--compare='):]
    elif arg.startswith('--save='):
        save = arg[len('--save='):]
    elif arg.startswith('--'):
        raise ValueError('Unrecognized option: ' + arg)
    else:
        args.append(arg)

inner = (inner != 0)
if inner:
    offset = -2127.
else:
    offset = 0.

plot_dir = plot
make_plots = (plot != '')

image_dir = images
make_images = (images != '')

if type(axis) == int:
    axes = [axis]
else:
    axes = axis

if len(args) < 3:
    raise IOError('missing arg(s)')

# Parameters of radial profiles
SMOOTH = 1000.  # km

validate = (validate != 0)
profile_flag = (profile != 0)

if profile_flag:
    a = np.fromfile(profile, sep=' ')
    a = a.reshape(a.size/2,2)
    profile_x = a[:,0]
    profile_y = a[:,1]

    validation_if = Tabulation(profile_x, profile_y)

if compare:
    a = np.fromfile(compare, sep=' ')
    a = a.reshape(a.size/2,2)
    compare_x = a[:,0]
    compare_y = a[:,1]

################################################################################
# Initialization
################################################################################

S_STEP = 10.**(1./3.)
RHO_VALUES = [1./S_STEP, 1., S_STEP]
Q_VALUES = [2.5, 3.0, 3.5]
BODY_TYPES = ['nominal', 'charon', 'kbo']
BODY_NAMES = {905:'styx', 902:'nix', 904:'kerberos', 903:'hydra'}

f = open(production)
production_0 = pickle.load(f)
f.close()
# indexing is [case][speed][rho][q]

# Identify body from first arg
body_id = int(args[0])

if body_id > 900:
    body_name = BODY_NAMES[body_id]
else:
    body_name = prefix + '-' + ('%04d' % body_id)

# Get destination dir from second arg
dest = args[1]

################################################################################
# Load the grids
################################################################################

big_dict = load_tree(args[2:], bodies=[body_id], verbose=True, drag=drag,
                               match=pattern)

big_dict = big_dict[body_id]
# indexing is [body_type][case][speed][beta]

# Insert the missing cases for an inner moon
if inner:
    for body_type in big_dict.keys():
        case_dict = big_dict[body_type]
        for case in (2,3,4):
            if case not in case_dict:
                case_dict[case] = case_dict[1]

################################################################################
# Scale each set of six size bins to normal I/F
################################################################################

if validate:
    bmin = -6
    bmax =  3
else:
    bmin = -4
    bmax =  1

save_bmin = -4
save_bmax =  1

for case in (1,2,3,4):
  for speed in (1,2):
    for body_type in big_dict.keys():

      # Skip combinations that have not been loaded
      try:
        beta_dicts = big_dict[body_type][case][speed]
      except KeyError:
        continue

      first_dict = beta_dicts.values()[0]

      shape = first_dict['grid'].shape
      xcell = first_dict['xsize']
      ycell = first_dict['ysize']
      zcell = first_dict['zsize']
      albedo = first_dict['albedo']

      info = np.array(shape + (xcell, ycell, zcell))

      for rho in RHO_VALUES:
        irho = int(rho)

        headers_with_grids = reindex_betas(beta_dicts, rho, bmin, bmax)

        for q in Q_VALUES:
            filename = '%s_speed%d_q%3.1f_case%d_%s_rho%d' % (
                          body_name, speed, q, case, BODY_TYPES[body_type],
                          irho)
            filepath = os.path.join(dest, filename + '.npz')
            if os.path.exists(filepath): continue

            # Calibrate the 3-D grids to number per km**3
            density = density_3d(headers_with_grids, q,
                                 production_0[case][speed][rho][q])

            # Convert the dictionary to a 3-D array
            merged = np.empty((save_bmax - save_bmin + 1,) + density[0].shape,
                              dtype='|f4')
            for b in range(save_bmin, save_bmax+1):
                merged[b - save_bmin] = density[b]

            # Write numpy file
            np.savez(filepath, info=info, density=merged)

            # Sum the image
            image = np.zeros(shape[:2])
            for b in density:
              scaling = np.pi * S_STEP**(2*b) * 1.e-12 * albedo * zcell
              image += scaling * density[b].sum(axis=-1)

            # Print reverse-calibrated info
            if validate:
              model = image_from_profile(validation_if, first_dict)
              mask = (model > 0.)

              if SMOOTH > 0.:
                image = gaussian_filter(image, (SMOOTH/xcell, SMOOTH/ycell))

              print image.max(), (image[mask]/model[mask]).max(), filename

            else:
              print image.max(), filename

            # Create radial profile if needed
            if make_plots or save:
                image = if_norm_2d(headers_with_grids, q,
                                   production_0[case][speed][rho][q])

                rmax = max(xcell * shape[0]/2., ycell * shape[1]/2.)
                dr = min(xcell, ycell) / 2.
                samples = int(np.ceil(rmax/dr))
                profile = profile_from_grid(image, first_dict, dr, samples,
                                            offset=offset)

            # Save plot if necessary
            if make_plots:
                pylab.close()

                if profile_flag:
                    pylab.plot(profile_x/1000., profile_y, 'r')

                if compare:
                    pylab.plot(compare_x/1000., compare_y, 'g')

                pylab.plot(dr/1000. * np.arange(samples), profile, 'k', lw=2)
                pylab.xlabel('Radius (1000 km)')
                pylab.ylabel('Normal I/F')

                filepath = os.path.join(plot_dir, filename) + '_profile.pdf'
                pylab.savefig(filepath)

            # Save images if necessary:
            if make_images:
              for axis in axes:
                image = if_norm_2d(headers_with_grids, q,
                                   production_0[case][speed][rho][q], axis=axis)
                filepath = os.path.join(image_dir, filename) + \
                                '_' + 'xyz'[axis] + 'image.png'

                save_image(image, filepath, 0., None)

            # Save a tabulation if necessary
            if save:
                f = open(os.path.join(save,filename) + '_profile.tab', 'w')

                for i in range(samples):
                    f.write('%5.0f ' % (i * dr))
                    f.write('%10.4e\n' % profile[i])

                f.close()

################################################################################
