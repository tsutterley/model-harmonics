#!/usr/bin/env python
"""
racmo_smb_cumulative.py
Written by Tyler Sutterley (07/2026)
Reads RACMO datafiles to calculate cumulative anomalies in derived surface
    mass balance products

CALLING SEQUENCE:
    python racmo_smb_cumulative.py --product smb --verbose <path_to_racmo_file>

COMMAND LINE OPTIONS:
    -P X, --product X: RACMO SMB product to calculate
    --mean: Start and end year of mean
    -G, --gzip: netCDF4 file is locally gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Local permissions mode of the directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 07/2026: output using structured dictionary with netCDF4 parameters
    Updated 06/2025: generalize the regular expression for model region
        can remap netCDF4 variable names for specified product
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: deprecation fixes for regular expressions
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: complete rewrite of program
        dropped old RACMO ascii file read portions
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: concatenated numpy arrays as input to nearest neighbors
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: can use getopt to set the working directory
    Updated 01/2018: using scikit-learn nearest neighbors for quickly finding
        mapping between input model grid and output global grid
    Updated 10/2017: new way of selecting between model versions
    Updated 09/2017: added RACMO2.3 XANT27 1979-2016
    Updated 02/2017: base_dir as input to racmo functions
    Updated 07/2016: added RACMO2.3 1958-2015 model, using netCDF4-python
    Updated 06/2016: using __future__ print function
    Updated 12/2015: RACMO2.3 XANT27 run for 1979-2014
    Updated 06/2015: added new PRODUCTS including derived Rainfall product
        saving indices for converting from RACMO grid to global grid
    Updated 04/2015: added HDF5 output option
    Updated 03/2015: updated for new test_ant27 program
    Updated 06/2014: can run from command line and parameter file
    Updated 05/2014: code generalizations and updates
    Updated 02/2014: more general code updates
    Updated 01/2014: create date files. add new RACMO 2.3 model
    Written 10/2011
"""

from __future__ import print_function

import sys
import re
import copy
import gzip
import uuid
import time
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: read and cumulative RACMO SMB SMB estimates
def racmo_smb_cumulative(
    model_file, PRODUCT, RANGE=None, VARIABLE=None, GZIP=False, MODE=0o775
):
    # get logger
    logger = logging.getLogger(__name__)
    # RACMO SMB model file
    model_file = pathlib.Path(model_file).expanduser().absolute()
    # try to extract region and version from filename
    R1 = re.compile(r'[P]?[XF]?(ANT|GRN|PEN|ASE)(\d+)', re.VERBOSE)
    R2 = re.compile(r'(RACMO\d+(\.\d+)?(p\d+)?)', re.VERBOSE)
    REGION = R1.search(model_file.name).group(0)
    VERSION = R2.search(model_file.name).group(0)
    # RACMO products
    racmo_products = {}
    racmo_products['precip'] = 'Precipitation'
    racmo_products['rainfall'] = 'Rainfall'
    racmo_products['refreeze'] = 'Meltwater Refreeze'
    racmo_products['runoff'] = 'Meltwater Runoff'
    racmo_products['smb'] = 'Surface Mass Balance'
    racmo_products['sndiv'] = 'Snow Drift Erosion'
    racmo_products['snowfall'] = 'Snowfall'
    racmo_products['snowmelt'] = 'Snowmelt'
    racmo_products['subl'] = 'Sublimation'
    # check if variable is remapped
    if VARIABLE is None:
        VARIABLE = copy.copy(PRODUCT)

    # Open the RACMO SMB NetCDF file for reading
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(str(model_file), mode='r') as f:
            f_in = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # Open the RACMO NetCDF file for reading
        f_in = netCDF4.Dataset(model_file, mode='r')

    # Output NetCDF file information
    logger.info(str(model_file))
    logger.info(list(f_in.variables.keys()))

    # dictionary with input/output data
    fd = {}

    # dictionary defining output structure
    struct = dict(
        dimensions=('time', 'rlat', 'rlon'),
        variables={
            PRODUCT: ('time', 'rlat', 'rlon'),
            'lat': ('rlat', 'rlon'),
            'lon': ('rlat', 'rlon'),
            'rotated_pole': (),
        },
    )

    # dictionary defining file-level and variable attributes
    attrs = dict(ROOT={})
    # file-level attributes to retrieve
    root_attributes = [
        'comment',
        'Domain',
        'Experiment',
        'source',
        'title',
    ]
    # variable attributes to retrieve
    variable_attributes = [
        'axis',
        'calendar',
        'description',
        'grid_mapping',
        'long_name',
        'standard_name',
        'units',
    ]
    # global attributes of NetCDF file
    attrs['ROOT'].update(ncdf_attributes(f_in, root_attributes))
    # output custom file-level attributes
    attrs['ROOT']['description'] = (
        f'Cumulative anomalies in {VERSION} {REGION} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}'
    )
    attrs['ROOT']['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    attrs['ROOT']['lineage'] = model_file.name

    # get data and attributes for each netCDF variable
    for key in [VARIABLE, 'lon', 'lat', 'rlon', 'rlat', 'time']:
        # output variable name (check if remapped)
        var = PRODUCT if (key == VARIABLE) else key
        # remove singleton dimensions
        fd[var] = np.squeeze(f_in.variables[key][:].copy())
        # get applicable attributes for variable
        attrs[var] = ncdf_attributes(f_in.variables[key], variable_attributes)

    # get the fill value for the variable
    fill_value = f_in.variables[VARIABLE].getncattr('_FillValue')
    # copy variable and attributes for projection
    key = 'rotated_pole'
    fd[key] = np.byte()
    attrs[key] = ncdf_attributes(f_in.variables[key], f_in[key].ncattrs())

    # parse date string within netCDF4 file
    date_string = attrs['time']['units']
    epoch1, to_secs = gravtk.time.parse_date_string(date_string)
    # calculate Julian day by converting to MJD and adding offset
    JD = 2400000.5 + gravtk.time.convert_delta_time(
        to_secs * fd['time'],
        epoch1=epoch1,
        epoch2=(1858, 11, 17, 0, 0, 0),
        scale=1.0 / 86400.0,
    )
    # convert from Julian days to calendar dates
    YY, MM, DD, hh, mm, ss = gravtk.time.convert_julian(JD, FORMAT='tuple')
    # convert from calendar dates to year-decimal
    TIME = gravtk.time.convert_calendar_decimal(
        YY, MM, day=DD, hour=hh, minute=mm, second=ss
    )

    # copy data to masked array
    DATA = np.ma.masked_equal(fd[PRODUCT].copy(), fill_value)
    # input shape of RACMO data
    nt, ny, nx = np.shape(DATA)

    # calculate mean period for RACMO
    (tt,) = np.nonzero((TIME >= RANGE[0]) & (TIME < (RANGE[1] + 1)))
    # cumulative mass anomalies calculated by removing mean balance flux
    MEAN = np.mean(DATA.data[tt, :, :], axis=0)
    # indices of specified ice mask at the first slice
    i, j = np.nonzero(~DATA.mask[0, :, :])
    valid_count = np.count_nonzero(~DATA.mask[0, :, :])
    # allocate for output variable
    fd[PRODUCT] = np.ma.zeros((nt, ny, nx), fill_value=DATA.fill_value)
    fd[PRODUCT].mask = DATA.mask | np.isnan(DATA.data)
    CUMULATIVE = np.zeros((valid_count))
    # calculate output cumulative anomalies for variable
    for t in range(nt):
        # convert mass flux from yearly rate and
        # calculate cumulative anomalies at time t
        CUMULATIVE += DATA.data[t, i, j] - MEAN[i, j]
        fd[PRODUCT].data[t, i, j] = CUMULATIVE.copy()
    # replace masked values with fill value
    fd[PRODUCT].data[fd[PRODUCT].mask] = fd[PRODUCT].fill_value

    # Output NetCDF filename
    FILE = f'{VERSION}_{REGION}_{PRODUCT.upper()}_cumul.nc'
    output_file = model_file.with_name(FILE)
    logger.info(str(output_file))

    # output RACMO data file with cumulative data
    if GZIP:
        # write data to in-memory netCDF4 file
        nc_buffer = mdlhmc.spatial.to_netCDF4(
            uuid.uuid4().hex,
            fd,
            attrs,
            struct,
            mode='w',
            memory=True,
        )
        # write RACMO data file to gzipped file
        with gzip.open(str(output_file), 'wb') as f:
            f.write(nc_buffer)
    else:
        # write data to netCDF4 file
        mdlhmc.spatial.to_netCDF4(
            output_file,
            fd,
            attrs,
            struct,
            mode='w',
        )
    # change the permissions mode
    output_file.chmod(mode=MODE)


# PURPOSE: get attributes for netCDF4 files and variable
def ncdf_attributes(nc, attributes_list):
    # get logger
    logger = logging.getLogger(__name__)
    # output dictionary of attributes
    attributes = {}
    # for each attribute to try to get
    for att_name in attributes_list:
        rx = re.compile(att_name, re.IGNORECASE)
        try:
            # use case-insensitive regex to find attribute name
            (ncattr,) = [s for s in nc.ncattrs() if rx.match(s)]
            att_val = nc.getncattr(ncattr)
        except Exception as exc:
            ncvar = getattr(nc, 'name', 'ROOT')
            logger.debug(f'Attribute {att_name} not found in {ncvar}')
            continue
        # strip whitespace for string attributes
        if isinstance(att_val, str):
            att_val = att_val.strip()
        # add attribute to dictionary
        attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates cumulative anomalies of RACMO
            surface mass balance products
            """
    )
    # command line parameters
    parser.add_argument(
        'infile', type=pathlib.Path, help='RACMO SMB file to run'
    )
    # products from SMB model
    choices = [
        'precip',
        'rainfall',
        'refreeze',
        'runoff',
        'smb',
        'sndiv',
        'snowfall',
        'snowmelt',
        'subl',
    ]
    parser.add_argument(
        '--product',
        '-P',
        type=str,
        metavar='PRODUCT',
        default='smb',
        choices=choices,
        help='RACMO SMB product to calculate',
    )
    parser.add_argument(
        '--variable',
        '-N',
        type=str,
        default=None,
        help='RACMO netCDF variable name for product',
    )
    # start and end years to run for mean
    parser.add_argument(
        '--mean',
        '-m',
        metavar=('START', 'END'),
        type=int,
        nargs=2,
        default=[1980, 1995],
        help='Start and end year range for mean',
    )
    # print information about each input and output file
    parser.add_argument(
        '--verbose',
        '-V',
        action='count',
        default=0,
        help='Verbose output of processing run',
    )
    # netCDF4 files are gzip compressed
    parser.add_argument(
        '--gzip',
        '-G',
        default=False,
        action='store_true',
        help='netCDF4 file is locally gzip compressed',
    )
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logger = gravtk.utilities.build_logger(
        __name__, level=loglevels[args.verbose]
    )

    # run program
    racmo_smb_cumulative(
        args.infile,
        args.product,
        VARIABLE=args.variable,
        RANGE=args.mean,
        GZIP=args.gzip,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
