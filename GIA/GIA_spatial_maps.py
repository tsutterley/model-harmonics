#!/usr/bin/env python
u"""
GIA_spatial_maps.py
Written by Tyler Sutterley (09/2023)
Calculates spatial maps of Glacial Isostatic Adjustment (GIA)

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    -l X, --lmax X: maximum spherical harmonic degree
    -R X, --radius X: Gaussian smoothing radius (km)
    -G X, --gia X: GIA model type to read
        IJ05-R2: Ivins R2 GIA Models
        W12a: Whitehouse GIA Models
        SM09: Simpson/Milne GIA Models
        ICE6G: ICE-6G GIA Models
        Wu10: Wu (2010) GIA Correction
        AW13-ICE6G: Geruo A ICE-6G GIA Models
        Caron: Caron JPL GIA Assimilation
        ICE6G-D: ICE-6G Version-D GIA Models
        ascii: reformatted GIA in ascii format
        netCDF4: reformatted GIA in netCDF4 format
        HDF5: reformatted GIA in HDF5 format
    --gia-file X: GIA file to read
    -U X, --units X: output units
        1: cm of water thickness
        2: mm of geoid height
        4: microGal gravitational perturbation
        5: mbar equivalent surface pressure
        6: cm viscoelastic crustal deformation [Wahr 2000]
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Pythonic interface to the HDF5 binary data format
        (http://www.h5py.org/)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

PROGRAM DEPENDENCIES:
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    harmonics.py: data class for working with spherical harmonics
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 09/2023: public release of GIA spatial program
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of gravity toolkit
    Updated 11/2022: use f-strings for formatting verbose output
    Updated 05/2022: use argparse descriptions within documentation
        use GIA reference and citation output from GIA read program
    Updated 04/2020: updates to reading load love numbers
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2018: using python3 compatible octal and input
    Updated 10/2016: minor clean up of output filename concatenation
    Updated 05-06/2016: using __future__ print function, added MMAX strings
    Updated 02/2016: use getopt parameters to set number of PROCESSES to run in
        parallel, whether or not to output a log file, added new help module
    Updated 11/2015: will output unique logs with parameters and output_file
    Updated 05/2015: added parameter MMAX for MMAX != LMAX
    Updated 03/2015: added error handling with traceback
    Updated 12/2014: complete update to program
    Updated 02/2014: general code updates to match other programs
    Written 05/2013
"""
from __future__ import print_function

import sys
import os
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
import collections
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: converts GIA spherial harmonics to spatial maps
def GIA_spatial_maps(LMAX,
    LMIN=1,
    MMAX=None,
    RAD=0,
    GIA=None,
    GIA_FILE=None,
    UNITS=None,
    DDEG=None,
    INTERVAL=None,
    BOUNDS=None,
    DATAFORM=None,
    OUTPUT_DIRECTORY=None,
    MODE=0o775):

    # output attributes for spatial files
    attributes = collections.OrderedDict()
    attributes['product_type'] = 'gravity_field'
    attributes['title'] = 'Glacial Isostatic Adjustment (GIA) Spatial Data'
    # list object of output files for file logs (full path)
    output_files = []
    # default output directory is the GIA directory
    if not OUTPUT_DIRECTORY:
        OUTPUT_DIRECTORY = GIA_FILE.parent

    # file information
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]

    # Calculating the Gaussian smoothing for radius RAD
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
        gw_str = f'_r{RAD:0.0f}km'
        attributes['smoothing_radius'] = f'{RAD:0.0f} km'
    else:
        # else = 1
        wt = np.ones((LMAX+1))
        gw_str = ''

    # flag for spherical harmonic order
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    order_str = f'M{MMAX:d}' if (MMAX != LMAX) else ''
    # add attributes for LMAX and MMAX
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX

    # read arrays of kl, hl, and ll Love Numbers
    # these Love numbers are not used in the spatial calculation
    LOVE = gravtk.load_love_numbers(LMAX,
        LOVE_NUMBERS=0, REFERENCE='CF', FORMAT='class')
    # do not include the elastic component in the unit coefficients
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE, include_elastic=False)

    # input GIA spherical harmonic datafiles
    GIA_Ylms_rate = gravtk.gia(lmax=LMAX).from_GIA(GIA_FILE, GIA=GIA)
    FILE_PREFIX = f'{GIA_Ylms_rate.title}_'
    attributes['GIA'] = (str(GIA_Ylms_rate.citation), GIA_FILE.name)

    # Output spatial data object
    grid = gravtk.spatial()
    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
        grid.lon = -180 + dlon*np.arange(0,nlon)
        grid.lat = 90.0 - dlat*np.arange(0,nlat)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        grid.lon = np.arange(-180+dlon/2.0,180+dlon/2.0,dlon)
        grid.lat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        grid.lon = np.arange(minlon+dlon/2.0, maxlon+dlon/2.0, dlon)
        grid.lat = np.arange(maxlat-dlat/2.0, minlat-dlat/2.0, -dlat)
        nlon = len(grid.lon)
        nlat = len(grid.lat)

    # Computing plms for converting to spatial domain
    theta = (90.0-grid.lat)*np.pi/180.0
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # output spatial units
    # dfactor computes the degree dependent coefficients
    # do not include the elastic component in the unit coefficients
    factors = gravtk.units(lmax=LMAX).harmonic(*LOVE, include_elastic=False)
    # 1: cmwe, centimeters water equivalent
    # 2: mmGH, millimeters geoid height
    # 4: micGal, microGal gravity perturbations
    # 5: mbar, millibars equivalent surface pressure
    # 6: cmCU, cm crustal uplift for viscoelastic uplift
    units = gravtk.units.bycode(UNITS)
    dfactor = factors.get(units)
    # output spatial units and descriptive units longname
    units_name, units_longname = gravtk.units.get_attributes(units)
    # add attributes for earth parameters
    attributes['earth_radius'] = f'{factors.rad_e:0.3f} cm'
    attributes['earth_density'] = f'{factors.rho_e:0.3f} g/cm^3'
    attributes['earth_gravity_constant'] = f'{factors.GM:0.3f} cm^3/s^2'
    # add attributes to output spatial object
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'
    grid.attributes['ROOT'] = attributes

    # converting harmonics to truncated, smoothed coefficients in units
    # combining harmonics to calculate output spatial fields
    Ylms = GIA_Ylms_rate.copy().convolve(dfactor*wt)
    # convert spherical harmonics to output spatial grid
    grid.data = gravtk.harmonic_summation(Ylms.clm, Ylms.slm,
        grid.lon, grid.lat, LMIN=LMIN, LMAX=LMAX,
        MMAX=MMAX, PLM=PLM).T
    grid.mask = np.zeros_like(grid.data, dtype=bool)

    # output files to ascii, netCDF4 or HDF5
    FILE = f'{FILE_PREFIX}{units}_L{LMAX:d}{order_str}{gw_str}.{suffix}'
    OUTPUT_FILE = OUTPUT_DIRECTORY.joinpath(FILE)
    grid.to_file(OUTPUT_FILE, format=DATAFORM, date=False,
        units=f'{units_name} yr^1', longname=units_longname)
    # set the permissions mode of the output files
    OUTPUT_FILE.chmod(mode=MODE)
    # add file to list
    output_files.append(OUTPUT_FILE)

    # return the list of output files
    return output_files

# PURPOSE: print a file log for the GIA spatial conversion
def output_log_file(input_arguments, output_files):
    # format: GIA_spatial_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GIA_spatial_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print output files
    logging.info('\n\nOUTPUT FILES:')
    for f in output_files:
        logging.info(f)
    # close the log file
    fid.close()

# PURPOSE: print a error file log for the GIA spatial conversion
def output_error_log_file(input_arguments):
    # format: GIA_spatial_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'GIA_spatial_failed_run_{0}_PID-{1:d}.log'.format(*args)
    # create a unique log and open the log file
    DIRECTORY = pathlib.Path(input_arguments.output_directory)
    fid = gravtk.utilities.create_unique_file(DIRECTORY.joinpath(LOGFILE))
    logging.basicConfig(stream=fid, level=logging.INFO)
    # print argument values sorted alphabetically
    logging.info('ARGUMENTS:')
    for arg, value in sorted(vars(input_arguments).items()):
        logging.info(f'{arg}: {value}')
    # print traceback error
    logging.info('\n\nTRACEBACK ERROR:')
    traceback.print_exc(file=fid)
    # close the log file
    fid.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates spatial maps of Glacial
            Isostatic Adjustment (GIA)
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        gravtk.utilities.convert_arg_line_to_args
    # minimum spherical harmonic degree
    parser.add_argument('--lmin',
        type=int, default=1,
        help='Minimum spherical harmonic degree')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # Gaussian smoothing radius (km)
    parser.add_argument('--radius','-R',
        type=float, default=0,
        help='Gaussian smoothing radius (km)')
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['Caron'] = 'Caron JPL GIA Assimilation'
    models['ICE6G-D'] = 'ICE-6G Version-D GIA Models'
    models['ascii'] = 'reformatted GIA in ascii format'
    models['netCDF4'] = 'reformatted GIA in netCDF4 format'
    models['HDF5'] = 'reformatted GIA in HDF5 format'
    # GIA model type
    parser.add_argument('--gia','-G',
        type=str, metavar='GIA', choices=models.keys(),
        help='GIA model type to read')
    # full path to GIA file
    parser.add_argument('--gia-file',
        type=pathlib.Path,
        help='GIA file to read')
    # output units
    parser.add_argument('--units','-U',
        type=int, default=1, choices=[1,2,4,5,6],
        help='Output units')
    # output grid parameters
    parser.add_argument('--spacing',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval',
        type=int, default=2, choices=[1,2,3],
        help=('Output grid interval '
            '(1: global, 2: centered global, 3: non-global)'))
    parser.add_argument('--bounds',
        type=float, nargs=4, metavar=('lon_min','lon_max','lat_min','lat_max'),
        help='Bounding box for non-global grid')
    # input data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input/output data format')
    parser.add_argument('--output-directory','-O',
        type=pathlib.Path,
        help='Output directory for spatial files')
    # Output log file for each job in forms
    # GIA_spatial_run_2002-04-01_PID-00000.log
    # GIA_spatial_failed_run_2002-04-01_PID-00000.log
    parser.add_argument('--log',
        default=False, action='store_true',
        help='Output log file for each job')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # try to run the analysis with listed parameters
    try:
        info(args)
        # run algorithm with parameters
        output_files = GIA_spatial_maps(args.lmax,
            LMIN=args.lmin,
            MMAX=args.mmax,
            RAD=args.radius,
            GIA=args.gia,
            GIA_FILE=args.gia_file,
            UNITS=args.units,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            DATAFORM=args.format,
            OUTPUT_DIRECTORY=args.output_directory,
            MODE=args.mode)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
        if args.log:# write failed job completion log file
            output_error_log_file(args)
    else:
        if args.log:# write successful job completion log file
            output_log_file(args, output_files)

# run main program
if __name__ == '__main__':
    main()
