################################################################################
# nhdust.py
#
# Tools for manipulating dust models. Previously headers.py.
#
# version 1.0, 11/20/2014, Mark Showalter
# version 1.1, 12/8/2014, Small revision to formatting of beta in file names
# version 1.2, 12/23/2014, Allow floats to appear as ints in the header file;
#   renamed nhdust.py, added tools for grid manipulations including calibration.
# version 1.3, 2/17/2015, Added capacity to deal with drag.
# version 1.4, 4/6/2015, Added graphing features for individual integrations.
################################################################################

import os
import numpy as np
import pylab
import re
import glob
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import zoom as zoom_image

try:
    from PIL import Image
except ImportError:
    pass

################################################################################
# Global constants
################################################################################

HEADER_NAMES_TYPES = [
    ('imax','i'),
    ('jmax','i'),
    ('kmax','i'),
    ('xsize','f'),
    ('ysize','f'),
    ('zsize','f'),
    ('body_id','i'),
    ('body_type','i'),
    ('mass','f'),
    ('radius','f'),
    ('density','f'),
    ('albedo','f'),
    ('v_esc','f'),
    ('case','i'),
    ('speed','i'),
    ('beta','f'),
    ('duration','f'),
    ('time_step','f'),
    ('ejected','i'),
    ('retained','i'),
    ('drag', 'i')       # optional, for backward compatibility
]

S_STEP = 10.**(1./3.)
LOG_S_STEP = np.log(S_STEP)

BODY_NAMES = ['styx', 'nix', 'kerberos', 'hydra']
BODY_IDS = [905, 902, 904, 903]
BODY_NAME_DICT = {905: 'styx', 902: 'nix', 904: 'kerberos', 903: 'hydra'}
BODY_TYPES = ['nominal', 'kbo', 'charon']   # 0, 1, 2
RHO_VALUES = [1./S_STEP, 1., S_STEP]  # g/cm**3
Q_VALUES = [2.5, 3.0, 3.5]

SMOOTH = 1000.      # For Gaussian smoothing of 2-D image during calibration

PLUTO_XCENTER = -2127.

################################################################################
# Header manipulation routines
################################################################################

def validate_dict(header_dict):
    """Validates a dictionary of parameters or raises AssertionError."""

    # Deal with drag parameter
    if 'drag' not in header_dict:
        ldict = len(HEADER_NAMES_TYPES) - 1
    else:
        ldict = len(HEADER_NAMES_TYPES)

    # Check types
    for (name, typechar) in HEADER_NAMES_TYPES[:ldict]:
        value = header_dict[name]
        if typechar == 'i':
            if not isinstance(value, int):
                raise ValueError('Parameter %s should be of type int: %s' %
                                 (name, type(value).__name__))
            header_dict[name] = int(value)

        elif typechar == 'f':
            if not isinstance(value, (float,int)):
                raise ValueError('Parameter %s should be of type float: %s' %
                                 (name, type(value).__name__))
            header_dict[name] = float(value)

    assert len(header_dict) == ldict    # no extra parameters

    # Continue validation
    body_id = header_dict['body_id']
    assert body_id in range(902,906) + range(20)

    if body_id > 900:
        assert header_dict['body_type'] == 0
    else:
        assert header_dict['body_type'] in (1,2)

    assert header_dict['case'] in (1,2,3,4)

    assert header_dict['speed'] in (1,2)

    if 'drag' in header_dict:
        assert header_dict['drag'] in range(10)

def read_header(filepath):
    """Read the parameters from a header file and return a dictionary."""

    # Read the file records
    f = open(filepath, 'r')
    recs = f.readlines()
    f.close()

    # Initialize the dictionary
    header_dict = {}

    # Read and validate the records
    for i in range(len(recs)):
        pair = recs[i].split('#')
        assert len(pair) == 2

        name = pair[1].strip()
        value = eval(pair[0].strip())
        header_dict[name] = value

    # Fill in drag from the file path if necessary
    if 'drag' not in header_dict:
        k = filepath.find('drag')
        if k >= 0:
            header_dict['drag'] = int(filepath[k+4])
        else:
            header_dict['drag'] = 1     # default value

    validate_dict(header_dict)

    # Return the dictionary
    return header_dict

def write_header(filepath, header_dict):
    """Write a dictionary of parameters into a header file."""

    validate_dict(header_dict)

    # Deal with drag factor
    if 'drag' not in header_dict:
        ldict = len(HEADER_NAMES_TYPES) - 1
    else:
        ldict = len(HEADER_NAMES_TYPES)

    # Write the parameter file
    f = open(filepath, 'w')

    for (name, typechar) in HEADER_NAMES_TYPES[:ldict]:
        value = str(header_dict[name])
        f.write('%-11s # %s\n' % (value, name))

    f.close()

def file_prefix(header_dict, more='', body_prefix='new'):
    """Construct a file name from a dictionary, excluding the extension."""

    validate_dict(header_dict)

    fields = []

    # Field 1: body name
    body_id = header_dict['body_id']
    if body_id > 900:
        field = BODY_NAME_DICT[body_id]
    else:
        field = body_prefix + str(body_id)

    fields.append(field)

    # Field 2: body type
    body_type = header_dict['body_type']
    field = BODY_TYPES[body_type]
    fields.append(field)

    # Field 3: case
    field = 'case' + str(header_dict['case'])
    fields.append(field)

    # Field 4: speed distribution
    field = 'speed' + str(header_dict['speed'])
    fields.append(field)

    # Field 5: beta
    field = 'beta_%7.1e' % header_dict['beta']
    fields.append(field)

    # Field 6: drag (if not inside given suffix)
    if 'drag' not in more:
        if 'drag' in header_dict:
            field = 'drag' + str(header_dict['drag'])
        else:
            field = 'drag1'

        fields.append(field)

    # Additional fields
    if more:
        fields.append(more)

    return '_'.join(fields)

################################################################################
# Grid routines
################################################################################

def load_header_and_grid(filepath, alt_names=['gridC.array', 'gridC.array2',
                                              'model.array', 'model.array2'],
                         verbose=False):
    """Return a header dictionary and a grid given the full path to the header
    file.

    Input:
        filepath        complete path to the file.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             (header dictionary, associated 3-D grid)
    """

    header_dict = read_header(filepath)
    grid = load_grid(filepath, header_dict, alt_names, verbose)
    return (header_dict, grid)

def load_header_with_grid(filepath, alt_names=['gridC.array', 'gridC.array2',
                                               'model.array', 'model.array2'],
                          verbose=False):
    """Return a header dictionary with the associated 3-D grid embedded under
    the key 'grid'.

    Input:
        filepath        complete path to the file.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             (header dictionary, associated 3-D grid)
    """

    (header_dict, grid) = load_header_and_grid(filepath, alt_names, verbose)
    header_dict['grid'] = grid
    return header_dict

def load_grid(filepath, header_dict=None,
                        alt_names=['gridC.array', 'gridC.array2',
                                   'model.array', 'model.array2'],
                        verbose=False):
    """Return a 3-D grid given the full path to the header file.

    Input:
        filepath        complete path to the file.
        header          header dictionary associated with the file. If None,
                        the header is read first.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             the associated 3-D grid.
    """

    if verbose: print filepath

    if header_dict is None:
        header_dict = read_header(filepath)

    shape = (header_dict['imax'], header_dict['jmax'], header_dict['kmax'])

    # Break down the file path
    (prefix, ext) = os.path.splitext(filepath)

    # Look for the array as a binary dump
    arrayfile = prefix + '.array'
    if os.path.exists(arrayfile):
        return load_grid_type1(arrayfile, shape)

    # Look for the array as a sequence of indices and values
    arrayfile = prefix + '.array2'
    if os.path.exists(arrayfile):
        return load_grid_type2(arrayfile, shape)

    # Try the alternative filenames
    dir = os.path.split(filepath)[0]
    for alt in alt_names:
        arrayfile = os.path.join(dir, alt)
        if os.path.exists(arrayfile):
            if alt.endswith('2'):
                return load_grid_type2(arrayfile, shape)
            else:
                return load_grid_type1(arrayfile, shape)

    raise IOError('array file not found for ' + filepath)

def load_grid_type1(filepath, shape):
    """Load a grid from a binary dump file. IOError on failure."""

    return np.fromfile(filepath, dtype='<f4').reshape(shape)

def load_grid_type2(filepath, shape):
    """Load a grid from a tabulation by cell index. IOError on failure."""

    dtype = np.dtype([('i','<u2'),
                      ('j','<u2'),
                      ('k','<u2'),
                      ('x','<f4')])

    table = np.fromfile(filepath, dtype=dtype)

    grid = np.zeros(shape, dtype='f4')
    i = table['i'] - 1
    j = table['j'] - 1
    k = table['k'] - 1
    x = table['x']

    grid[i,j,k] = x

    return grid

################################################################################
# Subset routines
################################################################################

def load_header_with_grids(pattern, alt_names=['gridC.array', 'gridC.array2',
                                               'model.array', 'model.array2'],
                           verbose=False):
    """Return a list of dictionaries from header files that match a given match
    pattern.

    Input:
        pattern         a regular expression that matches all the subsets to be
                        combined.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             (header dictionary, associated 3-D grid)
    """

    header_dicts = []
    for filepath in glob.glob(pattern):
        try:
            header_dict = load_header_with_grid(filepath, alt_names, verbose)

        except IOError:
            print '**** missing grid file for ' + filepath
            continue

        except IndexError as e:
            print '**** index error in grid file for ' + filepath
            print e
            continue

        header_dicts.append(header_dict)

    return header_dicts

def load_subsets(filepath, alt_names=['gridC.array', 'gridC.array2',
                                      'model.array', 'model.array2'],
                 verbose=False):
    """Return a list of dictionaries for all the subsets that match at given
    filepath for one subset.

    Input:
        filepath        a path to a header file.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             (header dictionary, associated 3-D grid)
    """

    k = -6
    while (True):
        k = filepath.find('subset', k+6)
        if k < 0: break

        filepath = filepath[:k+6] + '*' + filepath[k+8:]

    return load_header_with_grids(filepath, alt_names, verbose)

def merge_subsets(filepath, alt_names=['gridC.array', 'gridC.array2',
                                      'model.array', 'model.array2'],
                  verbose=False):
    """Return a single merged dictionary for all the subsets that match at given
    filepath for one subset.

    Input:
        filepath        a path to a header file.
        alt_names       a list of one or more file names to try if the default
                        name is not found.
        verbose         True to list the path to each header file as it is
                        loaded.
    Return:             (header dictionary, associated 3-D grid)
    """

    header_dicts = load_subsets(filepath, alt_names, verbose)

    header_dict = header_dicts[0]
    grid = header_dict['grid']
    del header_dict['grid']

    for temp_dict in header_dicts[1:]:
        temp_grid = temp_dict['grid']
        del temp_dict['grid']
        assert temp_dict == header_dict

        grid += temp_grid

    grid /= len(header_dicts)
    header_dict['grid'] = grid

    return header_dict

################################################################################
# Model loader
################################################################################

def load_tree(paths, bodies=None, alt_names=['gridC.array', 'gridC.array2',
                                             'model.array', 'model.array2'],
                                  verbose=False, drag=1, match=None):
    """Walk one or more directory trees and loads every header and grid that it
    finds.

    Input:
        paths       a path string or list of path strings to walk.
        bodies      a body ID or a list of body IDs. Only bodies in this list
                    will be loaded. Use None to load all bodies.
        alt_names   an optional list of array file names to try to load if the
                    file with the default name is missing.
        verbose     True to print the names of header files as they are loaded.
        drag        which drag index (0, 1, 2, ...) to calibrate.
        match       a regular expression that must fall within the header file
                    in order to be included.

    Return:         a dictionary of dictionaries of headers, indexed by
                    [body_id][body_type][case][speed][beta]. Each grid is saved
                    inside the dictionary, indexed as 'grid'. If multiple
                    subsets of the grid are found, they are added together.
    """

    # Convert args to lists if necessary
    if type(paths) == str: paths = [paths]
    if isinstance(bodies, int): bodies = [bodies]

    if match is not None:
        match = re.compile(match)

    # Walk...
    big_dict = {}   # indexed [body_id][body_type][case][speed][beta]
    for path in paths:
      for (root, dirs, files) in os.walk(path, followlinks=True):
        for file in files:
            if not file.endswith('.header'): continue

            filepath = os.path.join(root, file)

            if match is not None:
                test = match.search(filepath)
                if test is None: continue

            # Extract header info
            header_dict = read_header(filepath)
            if header_dict['drag'] != drag: continue

            if verbose: print filepath

            this_dict = big_dict

            # Handle body_id
            body_id = header_dict['body_id']
            if bodies is not None and body_id not in bodies: continue
            if body_id not in this_dict:
                this_dict[body_id] = {}

            this_dict = this_dict[body_id]

            # Handle body_type
            body_type = header_dict['body_type']
            if body_type not in this_dict:
                this_dict[body_type] = {}

            this_dict = this_dict[body_type]

            # Handle case
            case = header_dict['case']
            if case not in this_dict:
                this_dict[case] = {}

            this_dict = this_dict[case]

            # Handle speed
            speed = header_dict['speed']
            if speed not in this_dict:
                this_dict[speed] = {}

            this_dict = this_dict[speed]

            # Handle beta, fill in header and grid
            beta = round_beta(header_dict['beta'])      # round for simplicity

            try:
                grid = load_grid(filepath, header_dict, alt_names)

            except IOError:
                print '**** missing grid file for ' + filepath
                continue

            except IndexError as e:
                print '**** index error in grid file for ' + filepath
                print e
                continue

            if beta not in this_dict:                   # if beta is new, fill
                header_dict['grid'] = grid
                this_dict[beta] = header_dict

            else:                                       # otherwise, merge
                # Confirm that header and grid match
                old_dict = this_dict[beta].copy()
                new_dict = header_dict.copy()
                del old_dict['ejected'], old_dict['retained'], old_dict['grid']
                del new_dict['ejected'], new_dict['retained']
                if old_dict != new_dict:
                    raise ValueError('header mismatch: ' + filepath)

                # Merge
                this_dict[beta]['ejected'] += header_dict['ejected']
                this_dict[beta]['retained'] += header_dict['retained']
                this_dict[beta]['grid'] += grid

    return big_dict

################################################################################
# Beta manipulations
################################################################################

def beta_value(rho, albedo, bin=1):
    """Return the value of beta given density, albedo and size bin.

    Input:
        rho             density in g/cm**3.
        albedo          geometric albedo
        bin             size bin. The size in mm is S_SSTEP**bin.

    Return:             beta
    """

    q_pr = 1 + 1.36 * albedo
    s = S_STEP**bin
    beta = 5.7e-4 * q_pr / rho / s  # s in mm, rho in g/cm**3

    return beta

def round_beta(beta):
    """Returns beta rounded to two significant digits, for use in indexing."""

    return float('%7.1e' % beta)

def rho_x_radius(beta, albedo):
    """Return the value of particle radius times density in g/cm**2, given
    beta and albedo.
    """

    q_pr = 1 + 1.36 * albedo
    return 5.7e-5 * q_pr / beta

################################################################################
# Grid interpolation
################################################################################

def header_with_grid(headers_with_grids, beta):
    """Return a header dictionary containing a grid, produced if necessary via
    interpolation from a dictionary of headers keyed by betas.
    """

    # See if it is already in the dictionary of headers
    rounded = round_beta(beta)
    if rounded in headers_with_grids:
        return headers_with_grids[rounded]

    # See if something matches within 5%
    keys = headers_with_grids.keys()
    keys.sort()

    diffs = np.abs(np.array(keys) / beta - 1.)
    test = np.where(diffs < 0.05)

    if len(test[0]) > 0:
        key = keys[test[0][0]]
        return headers_with_grids[key]

    # Very small betas can fall outside the tabulated limits
    if rounded <= keys[0]:
        new_header = headers_with_grids[keys[0]].copy()
        new_header['beta'] = beta
        return new_header

    # Extrapolation to larger beta should be flagged
    if rounded >= keys[-1]:
        print('Extrapolation for beta = %7.1e: (%7.1e)' % (beta, keys[-1]))
        new_header = headers_with_grids[keys[-1]].copy()
        new_header['beta'] = beta
        return new_header

    # Find the bounding keys
    key_array = np.array(keys)
    key_below = key_array[key_array < beta][-1]
    key_above = key_array[key_array > beta][0]

    # Interpolation should be flagged
    print('Interpolation for beta = %7.1e: (%7.1e, %7.1e)' %
                  (beta, key_below, key_above))

    # Interpolate and return
    header_below = headers_with_grids[key_below]
    header_above = headers_with_grids[key_above]

    beta_below = header_below['beta']
    beta_above = header_above['beta']

    frac = (beta - beta_below) / (beta_above - beta_below)
    new_header = header_below.copy()
    new_header['beta'] = beta
    new_header['grid'] =      frac  * header_above['grid'] + \
                         (1 - frac) * header_below['grid']

    return new_header

def reindex_betas(headers_with_grids, rho, bin_min=-4, bin_max=1):
    """Given a dictionary keyed by beta, return one keyed by size bin index.

    Note: This procedure does not confirm that all the grids match in shape and
    size.

    Input:
        headers_with_grids  a dictionary of header information with included
                            grids, keyed by beta.
        rho                 grain density in g/cm**3.
        bin_min             the smallest size bin required.
        bin_max             the largest size bin required.

    Return:                 a new dictionary with a complete set of bin indices
                            from bin_min to bin_max, inclusive.
    """

    albedo = headers_with_grids.values()[0]['albedo']

    new_betas = {}
    for b in range(bin_min, bin_max+1):
        beta = beta_value(rho, albedo, b)
        new_betas[b] = header_with_grid(headers_with_grids, beta)

    return new_betas

################################################################################
# Scaled models
################################################################################

def if_norm_2d(headers_with_grids, q, production_0=1., axis=-1):
    """Return a 2-D image which, when multiplied by production_0, is in units of
    normal I/F.

    Input:
        headers_with_grids  a dictionary of headers indexed by the size bin
                            index b. The 3-D grids are embedded within the
                            dictionaries, indexed by 'grid'.
        q                   the power-law index of the size distribution, 2.5,
                            3, or 3.5.
        production_0        the global surface production rate (number per km**3
                            per sec) for size bin p=0.
        axis                the axis to sum over: -1 for a top-down view;
                            0 or 1 for side views.

    Return:         a 2-D grid which, when multiplied by parameter production_0,
                    will be in units of normal I/F.
    """

    # Extract required header parameters
    header = headers_with_grids.values()[0]     # Needed parameters are the same
                                                # in every header
    n_esc = header['ejected']
    n_kept = header['retained']
    n_gen = n_esc + n_kept

    dx = header['xsize']
    dy = header['ysize']
    dz = header['zsize']
    assert (dx == dy) and (dy == dz)

    radius = header['radius']
    albedo = header['albedo']

    t_step = header['time_step']

    # Sum IF_norm / production_0 over each size bin and over the k-axis
    shape = (header['imax'], header['jmax'])
    if_norm = np.zeros(shape)

    scaling = ((production_0 * albedo * np.pi * 4.e-12 * np.pi * radius**2 *
                t_step) / (n_gen * dx * dy))

    for (b,header) in headers_with_grids.iteritems():
        grid = header['grid']
        if_norm += S_STEP**(b * (3. - q)) * grid.sum(axis=axis)

    if_norm *= scaling
    return if_norm

def density_3d(headers_with_grids, q, production_0=1.):
    """Return a dictionary of grids calibrated in units of number per km**3.

    Input:
        headers_with_grids  a dictionary of headers indexed by the size bin
                            index b. The 3-D grids are embedded within the
                            dictionaries, indexed as 'grid'.
        q                   the power-law index of the size distribution, 2.5,
                            3, or 3.5.
        production_0        the global surface production rate (number per km**3
                            per sec) for size bin p=0.

    Return:             a dictionary of 3-D grids, calibrated in units of number
                        per km**3, and indexed by size bin index b.
    """

    # Extract required header parameters
    header = headers_with_grids.values()[0]     # Needed parameters are the same
                                                # in every header
    n_esc = header['ejected']
    n_kept = header['retained']
    n_gen = n_esc + n_kept

    dx = header['xsize']
    dy = header['ysize']
    dz = header['zsize']

    radius = header['radius']

    t_step = header['time_step']

    # Scale the grids into units of km**(-3)
    scaling = ((production_0 * 4. * np.pi * radius**2 * t_step) /
               (n_gen * dx * dy * dz))
    density = {}
    for b in headers_with_grids:
        density[b] = scaling * S_STEP**(b * (1.-q)) * \
                     headers_with_grids[b]['grid']

    # Pile larger material into the upper bin
    b_max = np.max(headers_with_grids.keys())
    density[b_max] /= (1. - S_STEP**(1.-q))

    return density

################################################################################
# Radial profile operations
################################################################################

def profile_from_grid(grid, header_dict, dr, samples, smooth=0., offset=0.):
    """Create a radial profile from a 2-D or 3-D grid.

    Input:
        grid        2-D or 3-D array. If 3-D, a sum over the last axis is
                    performed.
        header_dict header dictionary describing the grid.
        dr          radial step to use in profile, in km.
        samples     number of samples in profile. The outer radial limit is
                    dr * samples.
        smooth      the amount by which to smooth the 2-D image before sampling,
                    in units of pixels. This is used as the standard deviation
                    of a Gaussian filter applied to the image. It can be a tuple
                    to support different smoothing along the x and y axes.
        offset      the offset of the ring center along the x-axis. Default 0.
                    Use -2127. for a Pluto-centered profile.
    Return:         a 1-D array containing the radial profile, after averaging
                    in longitude.
    """

    # Convert to 2-D
    if len(grid.shape) > 2:
        grid = grid.sum(axis=-1)

    # Smooth
    xcell = header_dict['xsize']
    ycell = header_dict['ysize']
    if smooth != 0.:
        grid = gaussian_filter(grid, (smooth/xcell, smooth/ycell))

    # Calculate extreme radii at the four corners of each cell
    imax = header_dict['imax']
    jmax = header_dict['jmax']
    x = xcell * (np.arange(imax+1) - imax/2.)[:,np.newaxis] - offset
    y = ycell * (np.arange(jmax+1) - jmax/2.)

    rstack = np.empty((4,) + (imax,jmax))
    rstack[0] = np.sqrt(x[:-1]**2 + y[:-1]**2)
    rstack[1] = np.sqrt(x[:-1]**2 + y[1:]**2)
    rstack[2] = np.sqrt(x[1:]**2  + y[:-1]**2)
    rstack[3] = np.sqrt(x[1:]**2  + y[1:]**2)

    rmin = rstack.min(axis=0)
    rmax = rstack.max(axis=0)

    # Find the mean among all the pixels that sample each selected radius
    profile = np.zeros(samples)
    for k in range(samples):
          r = k * dr
          mask = (r >= rmin) & (r <= rmax)
          if np.any(mask):
            profile[k] = np.mean(grid[mask])

    return profile

def image_from_profile(profile, header_dict):
    """Create a radial profile from a 2-D or 3-D grid.

    Input:
        profile     a tabulation of a ring property vs. radius.
        header_dict header dictionary describing the image to simulate.

    Return:         a 2-D array constructed by mapping the profile onto a 2-D
                    grid.
    """

    # Construct a 2-D array of radius
    xcell = header_dict['xsize']
    ycell = header_dict['ysize']
    imax = header_dict['imax']
    jmax = header_dict['jmax']

    x = xcell * (np.arange(imax) - imax/2. + 0.5)[:,np.newaxis]
    y = ycell * (np.arange(jmax) - jmax/2. + 0.5)
    r = np.sqrt(x**2 + y**2)

    # Return the 2-D image
    return profile(r)

################################################################################
# Lifetime calculations
################################################################################

def lifetime_per_km3(header_w_grid, ejected_only=False):
    """Returns a 3-D grid providing mean particle lifetimes in units of years
    per cubic km."""

    grid = header_w_grid['grid']

    n_esc = header_w_grid['ejected']
    n_kept = header_w_grid['retained']
    n_gen = n_esc + n_kept

    if ejected_only:
        n_scale = n_esc
    else:
        n_scale = n_gen

    t_step = header_w_grid['time_step'] / (86400. * 365.25)

    dx = header_w_grid['xsize']
    dy = header_w_grid['ysize']
    dz = header_w_grid['zsize']

    return grid * (t_step / (n_scale * dx * dy * dz))

def lifetime_per_km2(header_w_grid, axis=2, ejected_only=False):
    """Returns a 2-D grid providing particle lifetimes in units of years per
    square km."""

    grid_3d = lifetime_per_km3(header_w_grid, ejected_only=ejected_only)

    dx = header_w_grid['xsize']
    dy = header_w_grid['ysize']
    dz = header_w_grid['zsize']
    assert (dx == dy) and (dx == dz)

    return np.sum(grid_3d, axis=axis) * dx

def lifetime_per_km2_profile(header_w_grid,
                             dr=100., samples=1000, smooth=0., offset=0.,
                             ejected_only=False):
    """Returns a 1-D grid providing particle lifetimes in units of years per
    square km2, but averaged into a radial profile."""

    grid_2d = lifetime_per_km2(header_w_grid, axis=2, ejected_only=ejected_only)

    return profile_from_grid(grid_2d, header_w_grid,
                             dr, samples, smooth, offset)

def lifetime_per_km(header_w_grid, dr=100., samples=1000, smooth=0.,
                                   ejected_only=False):
    """Returns a 1-D radial profile providing particle lifetimes in units of
    years per km. The area under this curve should be the mean particle
    lifetime.

    Input:
        header_w_grid
        dr          radial step to use in profile, in km.
        samples     number of samples in profile. The outer radial limit is
                    dr * samples.
        smooth      the amount by which to smooth the 2-D image before sampling,
                    in units of pixels. This is used as the standard deviation
                    of a Gaussian filter applied to the image. It can be a tuple
                    to support different smoothing along the x and y axes.
    """

    profile = lifetime_per_km2_profile(header_w_grid, dr, samples, smooth,
                                       ejected_only=ejected_only)
    circumference = 2. * np.pi * dr * np.arange(samples)

    return profile * circumference

def mean_lifetime(header_w_grid, ejected_only=False):
    """Returns the mean lifetime in years among all the particles in the
    grid."""

    grid = header_w_grid['grid']

    n_esc = header_w_grid['ejected']
    n_kept = header_w_grid['retained']
    n_gen = n_esc + n_kept

    if ejected_only:
        n_scale = n_esc
    else:
        n_scale = n_gen

    t_step = header_w_grid['time_step'] / (86400. * 365.25)

    return np.sum(grid) * (t_step / n_scale)

################################################################################
# Graphical output
################################################################################

def save_image(image, filename, lo=None, hi=None, zoom=1.):
    """Save an image file (not necessarily png) of a 2-D array.

    Input:
        image       a 2-D array.
        fileanme    the name of the output file, which should end with the type,
                    e.g., '.png' or '.jpg'
        lo          the array value to map to black; if None, then the minimum
                    value in the array is used.
        hi          the array value to map to white; if None, then the maximum
                    value in the array is used.
        zoom        the scale factor by which to enlarge the image, default 1.
    """

    if lo is None:
        lo = image.min()

    if hi is None:
        hi = image.max()

    if zoom != 1:
        image = zoom_image(image, zoom)

    scaled = (image.swapaxes(0,1)[::-1].copy() - lo) / float(hi - lo)
    bytes = (256.*scaled).clip(0,255).astype("uint8")

    im = Image.fromstring("L", (bytes.shape[1], bytes.shape[0]), bytes)
    im.save(filename)

def save_plot(profiles, filename, rmin=None, rmax=None, ylabel='',
              styles=[], weights=[]):
    """Creates a line plot showing one or more a radial ring profiles.

    Input:
        profiles    a list of (x,y) tuples or lists, containing the x and y
                    values of each plot as NumPy arrays. The x array is in km.
        filename    name of the output file. The extension (.pdf, .jpg, etc.)
                    defines the type of the file.
        rmin        optional minimum radius value for the plot, km.
        rmax        optional maximum radius value for the plot, km.
        styles      an optional list of style types for the lines, in the same
                    order as the profiles. A None entry uses the default.
        weights     an optional list of line weights, in the same order as the
                    profiles.
    """

    pylab.close()

    # Determine the inner and outer limits if necessary
    if rmin is None or rmax is None:
        xmin = float(np.min(profiles[0][0]))
        xmax = float(np.max(profiles[0][0]))

        for i in range(1, len(profiles)):
            xmin = min(xmin, np.min(profiles[i][0]))
            xmax = max(xmax, np.max(profiles[i][0]))

        if rmin is None: rmin = xmin
        if rmax is None: rmax = xmax

    # Plot each curve, properly masked
    for i in range(len(profiles)):
        x = np.asfarray(profiles[i][0])
        y = np.asfarray(profiles[i][1])

        style = None
        if i < len(styles): style = styles[i]

        lweight = 1
        if i < len(weights): lweight = weights[i]

        mask = (x >= rmin) & (x <= rmax)

        if style is None:
            pylab.plot(x[mask]/1000., y[mask], lw=lweight)
        else:
            pylab.plot(x[mask]/1000., y[mask], style, lw=lweight)

        pylab.xlabel('Radius (1000 km)')
        if ylabel: pylab.ylabel(ylabel)

        pylab.savefig(filename)

################################################################################
