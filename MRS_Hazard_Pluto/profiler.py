#!/usr/bin/env python
################################################################################
# profiler.py
#
# Syntax:
#   python profiler.py --param=value ... directories...
#
# This program generates PDF diagrams of ring radial profiles for selected
# header/array files within a directory tree. Each plot shows curves for all of
# the available subsets.
#
# Inputs:
#   directories...      one or more directory trees that will be searched for
#                       header files and associated array and array2 files.
#
# Options:
#   --destination=str   directory path pointing to the destination of the output
#                       files.
#   --pattern=string    an optional shell match pattern to define which header
#                       files are used to generate ring profiles. This pattern
#                       should not specify a subset number. Example:
#                           --pattern=hydra_*drag0*
#   --verbose=n         1 to list the names of PDF files as they are generated;
#                       0 to supress this information. Default 1.
#   --subset=n          the number of one subset that exists. Default is 1. If
#                       only subset 0 exists, change to 0.
#   --inner=n           0 for a barycenter-centered ring; 1 for a Pluto-centered
#                       ring. Default is 0.
#   --images=dir        an optional directory in which to place PNG images
#                       derived from an average of all the subsets. Blank to
#                       suppress image generation. Default is blank.
#   --axis=n            the axis number over which the 3-D grid is collapsed to
#                       create the image. Default 2, for a sum over the z-axis
#                       producing a view of the x-y plane.
#
# Mark Showalter, mshowalter@seti.org
# Version 1.0: 4/6/15: Initial release.
# Version 1.1: 4/9/15: Added --inner option.
# Version 1.2: 4/14/15: Added image options.
################################################################################

import sys
import os
import numpy as np
import fnmatch
from nhdust import *

# Parameters of radial profiles
DR = 100.
SAMPLES = 1000

################################################################################
# Interpret command arguments
################################################################################

destination = './'
pattern = '*'
verbose = 1
subset = 1
inner = 0
images = ''
axis = 2

args = []

for arg in sys.argv[1:]:
    if arg.startswith('--destination='):
        destination = arg[len('--destination='):]
    elif arg.startswith('--pattern='):
        pattern = arg[len('--pattern='):]
    elif arg.startswith('--verbose='):
        exec(arg[2:])
    elif arg.startswith('--subset='):
        exec(arg[2:])
    elif arg.startswith('--inner='):
        exec(arg[2:])
    elif arg.startswith('--images='):
        images = arg[len('--images='):]
    elif arg.startswith('--axis='):
        exec(arg[2:])
    elif arg.startswith('--'):
        raise ValueError('Unrecognized option: ' + arg)
    else:
        args.append(arg)

if len(args) < 1:
    raise IOError('missing arg(s)')

verbose = (verbose != 0)
subset_match = 'subset%02d' % subset

if inner:
    offset = -2127.
else:
    offset = 0.

image_dir = images

################################################################################
# Generate profiles...
################################################################################

for path in args:
  for (root, dirs, files) in os.walk(path, followlinks=True):
    for file in files:

        (body, ext) = os.path.splitext(file)
        if ext != '.header': continue
        if subset_match not in body: continue
        if not fnmatch.fnmatch(body, pattern): continue
        filepath = os.path.join(root, file)

        subsets = load_subsets(filepath)
        header = subsets[0]

        dr = header['xsize'] / 2.
        samples = header['imax']
        x = dr * np.arange(samples)

        profiles = []
        images = []
        for header in subsets:
            y = lifetime_per_km2_profile(header, dr, samples, offset=offset)
            profiles.append((x,y))

            if image_dir:
                image = lifetime_per_km2(header, axis=axis)
                images.append(image)

        pdf_file = body[:-len('_subset01')] + '.pdf'
        outfile = os.path.join(destination, pdf_file)
        if verbose: print pdf_file

        save_plot(profiles, outfile, ylabel='Lifetime (years per square km)')

        if image_dir:
            png_file = pdf_file[:-4] + '.png'
            outfile = os.path.join(image_dir, png_file)
            save_image(np.mean(images, axis=0), outfile, 0.)

################################################################################
