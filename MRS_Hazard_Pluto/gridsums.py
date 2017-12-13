#!/usr/bin/env python
################################################################################
# gridsums.py
#
# Syntax:
#   python gridsums.py --param=value ... directories... output
#
# This program saves a dictionary mean grid lifetimes for selected header/array
# files within a directory tree. Each plot shows curves for all of the available
# subsets.
#
# Inputs:
#   directories...      one or more directory trees that will be searched for
#                       header files and associated array and array2 files.
#   output              name of the output directory, with extension .pickle
#
# Options:
#   --pattern=string    an optional shell match pattern to define which header
#                       files are used to generate ring profiles. This pattern
#                       should not specify a subset number. Example:
#                           --pattern=hydra_*drag0*
#   --subset=n          the number of one subset that exists. Default is 1. If
#                       only subset 0 exists, change to 0.
#
# Mark Showalter, mshowalter@seti.org
# Version 1.0: 4/12/15 Initial version.
################################################################################

import sys
import os
import numpy as np
import fnmatch
from nhdust import *
import pickle

################################################################################
# Interpret command arguments
################################################################################

destination = './'
pattern = '*'
verbose = 1
subset = 1
inner = 0

args = []

for arg in sys.argv[1:]:
    if arg.startswith('--pattern='):
        pattern = arg[len('--pattern='):]
    elif arg.startswith('--subset='):
        exec(arg[2:])
    elif arg.startswith('--'):
        raise ValueError('Unrecognized option: ' + arg)
    else:
        args.append(arg)

if len(args) < 2:
    raise IOError('missing arg(s)')

picklefile = args[-1]
assert picklefile.endswith('.pickle')

verbose = (verbose != 0)
subset_match = 'subset%02d' % subset

if inner:
    offset = -2127.
else:
    offset = 0.

################################################################################
# Generate profiles...
################################################################################

dictionary = {}

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

        lifetimes = []
        sums = []
        for header in subsets:
            y = mean_lifetime(header)
            lifetimes.append(y)
            sums.append(np.sum(header['grid']))

        info = body[:-len('_subset01')]
        print info, "%.5e" % np.mean(y), "%.5e" % np.mean(sums)
        dictionary[info] = (np.mean(y), np.mean(sums))

f = open(picklefile, 'wb')
pickle.dump(dictionary, f)
f.close()

################################################################################
