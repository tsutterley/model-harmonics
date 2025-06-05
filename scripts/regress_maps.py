#!/usr/bin/env python
u"""
regress_maps.py
Written by Tyler Sutterley (10/2024)

Reads in spatial files and fits a regression model at each grid point

INPUTS:
    path to input spatial file
    path to output spatial file

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output-directory X: output directory for spatial files
    -S X, --start X: starting GRACE month for time series regression
    -E X, --end X: ending GRACE month for time series regression
    -N X, --missing X: Missing GRACE/GRACE-FO months
    --spacing X: spatial resolution of output data (dlon,dlat)
    --interval X: output grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
        3: non-global grid (set with defined bounds)
    --bounds X: non-global grid bounding box (minlon,maxlon,minlat,maxlat)
    -F X, --format X: input/output data format
        ascii
        netCDF4
        HDF5
    --order X: regression fit polynomial order
    --cycles X: regression fit cyclical terms
    --log: Output log of files created for each job
    -V, --verbose: verbose output of processing run
    -M X, --mode X: Permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    time_series.regress.py: calculates trend coefficients using least-squares
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Written 10/2024
"""
from __future__ import print_function, division

import sys
import os
import time
import logging
import pathlib
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: keep track of threads
def info(args):
    logging.info(pathlib.Path(sys.argv[0]).name)
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# program module to run with specified parameters
def regress_maps(INPUT_FILE, OUTPUT_FILE,
        START=None,
        END=None,
        MISSING=None,
        DDEG=None,
        INTERVAL=None,
        BOUNDS=None,
        DATAFORM=None,
        ORDER=None,
        CYCLES=None,
        VERBOSE=0,
        MODE=0o775
    ):

    # operate on the input and output file paths
    INPUT_FILE = pathlib.Path(INPUT_FILE).expanduser().absolute()
    OUTPUT_FILE = pathlib.Path(OUTPUT_FILE).expanduser().absolute()
    # create output directory if currently non-existent
    OUTPUT_FILE.parent.mkdir(mode=MODE, parents=True, exist_ok=True)

    # Output Degree Spacing
    dlon,dlat = (DDEG[0],DDEG[0]) if (len(DDEG) == 1) else (DDEG[0],DDEG[1])
    # Output Degree Interval
    if (INTERVAL == 1):
        # (-180:180,90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):
        # (Degree spacing)/2
        nlon = np.int64(360.0/dlon)
        nlat = np.int64(180.0/dlat)
    elif (INTERVAL == 3):
        # non-global grid set with BOUNDS parameter
        minlon,maxlon,minlat,maxlat = BOUNDS.copy()
        lon = np.arange(minlon+dlon/2.0, maxlon+dlon/2.0, dlon)
        lat = np.arange(maxlat-dlat/2.0, minlat-dlat/2.0, -dlat)
        nlon = len(lon)
        nlat = len(lat)

    # read spatial file in data format
    if DATAFORM in ('ascii', 'netCDF4', 'HDF5'):
        # ascii (.txt)
        # netCDF4 (.nc)
        # HDF5 (.H5)
        dinput = gravtk.spatial().from_file(INPUT_FILE,
            format=DATAFORM, date=True, spacing=[dlon, dlat],
            nlat=nlat, nlon=nlon)
    elif DATAFORM in ('index-ascii', 'index-netCDF4', 'index-HDF5'):
        # read from index file
        _,dataform = DATAFORM.split('-')
        dinput = gravtk.spatial().from_index(INPUT_FILE,
            format=dataform, date=True, spacing=[dlon, dlat],
            nlat=nlat, nlon=nlon)
    # get the default attributes from the input data
    attributes = dinput.attributes['ROOT']
    units_name = dinput.attributes['data'].get('units')
    long_name = dinput.attributes['data'].get('long_name')

    # GRACE months to read
    if START is None:
        START = dinput.month[0]
    if END is None:
        END = dinput.month[-1]
    if MISSING is None:
        MISSING = []
    # create a list of months to fit
    months = sorted(set(np.arange(START,END+1)) - set(MISSING))
    # subset data to months
    dinput = dinput.subset(months)

    # Setting output parameters for each fit type
    coef_str = [f'x{o:d}' for o in range(ORDER+1)]
    if (ORDER == 0):# Mean
        fit_title = ['Mean']
    elif (ORDER == 1):# Trend
        fit_title = ['Constant','Trend']
    elif (ORDER == 2):# Quadratic
        fit_title = ['Constant','Linear','Quadratic']
    for i,c in enumerate(CYCLES):
        # check if fitting with semi-annual or annual terms
        if (c == 0.5):
            coef_str.extend(['SS','SC'])
            fit_title.extend(['Semi-Annual Sine', 'Semi-Annual Cosine'])
        elif (c == 1.0):
            coef_str.extend(['AS','AC'])
            fit_title.extend(['Annual Sine', 'Annual Cosine'])

    # Fitting seasonal components
    ncomp = (ORDER+1) + 2*len(CYCLES)
    # confidence interval for regression fit errors
    CONF = 0.95

    # Allocating memory for output variables
    out = dinput.zeros_like()
    out.data = np.zeros((nlat, nlon, ncomp))
    out.error = np.zeros((nlat, nlon, ncomp))
    out.mask = np.ones((nlat, nlon, ncomp),dtype=bool)
    out.time = np.arange(ncomp)

    # output attributes
    attr = dict(ROOT=attributes)
    attr['lon'] = {}
    attr['lon']['long_name'] = 'longitude'
    attr['lon']['units'] = 'degrees_east'
    attr['lat'] = {}
    attr['lat']['long_name'] = 'latitude'
    attr['lat']['units'] = 'degrees_north'
    attr['data'] = {}
    attr['data']['description'] = 'Model_fit'
    attr['data']['units'] = units_name
    attr['data']['long_name'] = long_name
    attr['data']['title'] = fit_title
    attr['error'] = {}
    attr['error']['description'] = 'Uncertainty_in_model_fit'
    attr['error']['units'] = units_name
    attr['error']['long_name'] = long_name
    attr['error']['confidence'] = 100*CONF
    attr['error']['title'] = fit_title
    attr['coefficients'] = dict(long_name='Regression_coefficients')
    # field mapping for output regression data
    field_mapping = {}
    field_mapping['lat'] = 'lat'
    field_mapping['lon'] = 'lon'
    field_mapping['data'] = 'data'
    field_mapping['error'] = 'error'
    field_mapping['time'] = 'coefficients'

    # Fit Significance
    # SSE: Sum of Squares Error
    # AIC: Akaike information criterion
    # BIC: Bayesian information criterion
    # R2Adj: Adjusted Coefficient of Determination
    for key in ['SSE','AIC','BIC','R2Adj']:
        setattr(out, key, np.zeros((nlat, nlon)))
        field_mapping[key] = key
    # output attributes for fit significance
    attr['SSE'] = dict(title='Sum of Squares Error')
    attr['AIC'] = dict(title='Akaike information criterion')
    attr['BIC'] = dict(title='Bayesian information criterion')
    attr['R2Adj'] = dict(title='Adjusted Coefficient of Determination')

    # calculate the regression coefficients and fit significance
    for i in range(nlat):
        for j in range(nlon):
            # Calculating the regression coefficients
            tsbeta = gravtk.time_series.regress(
                dinput.time, dinput.data[i,j,:],
                ORDER=ORDER, CYCLES=CYCLES, CONF=CONF)
            # save regression components
            for k in range(0, ncomp):
                out.data[i,j,k] = tsbeta['beta'][k]
                out.error[i,j,k] = tsbeta['error'][k]
                out.mask[i,j,k] = False
            # Fit significance terms
            # Degrees of Freedom
            nu = tsbeta['DOF']
            # Converting Mean Square Error to Sum of Squares Error
            out.SSE[i,j] = tsbeta['MSE']*nu
            out.AIC[i,j] = tsbeta['AIC']
            out.BIC[i,j] = tsbeta['BIC']
            out.R2Adj[i,j] = tsbeta['R2Adj']

    # output global attributes
    REFERENCE = f'Output from {pathlib.Path(sys.argv[0]).name}'
    # write to output file
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        out.to_ascii(OUTPUT_FILE, date=True, verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        out.to_netCDF4(OUTPUT_FILE, date=True, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attr,
            reference=REFERENCE)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        out.to_HDF5(OUTPUT_FILE, date=True, verbose=VERBOSE,
            field_mapping=field_mapping, attributes=attr,
            reference=REFERENCE)
    # change the permissions mode of the output file
    OUTPUT_FILE.chmod(mode=MODE)

# PURPOSE: print a file log for the regression
def output_log_file(input_arguments, output_files):
    # format: processing_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'processing_run_{0}_PID-{1:d}.log'.format(*args)
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

# PURPOSE: print a error file log for the regression
def output_error_log_file(input_arguments):
    # format: processing_failed_run_2002-04-01_PID-70335.log
    args = (time.strftime('%Y-%m-%d',time.localtime()), os.getpid())
    LOGFILE = 'processing_failed_run_{0}_PID-{1:d}.log'.format(*args)
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
        description="""Reads in spatial files and calculates the
            trends at each grid point following an input regression model
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # input and output file
    parser.add_argument('infiles',
        type=pathlib.Path, nargs='?',
        help='Input files')
    parser.add_argument('outfile',
        type=pathlib.Path, nargs='?',
        help='Output file')
    # start and end GRACE/GRACE-FO months
    parser.add_argument('--start','-S',
        type=int, default=None,
        help='Starting GRACE/GRACE-FO month for time series regression')
    parser.add_argument('--end','-E',
        type=int, default=None,
        help='Ending GRACE/GRACE-FO month for time series regression')
    parser.add_argument('--missing','-N',
        metavar='MISSING', type=int, nargs='+',
        help='Missing GRACE/GRACE-FO months')
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
    # input and output data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=choices,
        help='Input/output data format')
    # regression parameters
    # 0: mean
    # 1: trend
    # 2: acceleration
    parser.add_argument('--order',
        type=int, default=2,
        help='Regression fit polynomial order')
    # regression fit cyclical terms
    parser.add_argument('--cycles',
        type=float, default=[0.5,1.0], nargs='+',
        help='Regression fit cyclical terms')
    # Output log file for each job in forms
    # processing_run_2002-04-01_PID-00000.log
    # processing_failed_run_2002-04-01_PID-00000.log
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
        # run regress_maps algorithm with parameters
        output_files = regress_maps(args.infiles, args.outfile,
            START=args.start,
            END=args.end,
            MISSING=args.missing,
            DDEG=args.spacing,
            INTERVAL=args.interval,
            BOUNDS=args.bounds,
            DATAFORM=args.format,
            ORDER=args.order,
            CYCLES=args.cycles,
            VERBOSE=args.verbose,
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
            output_log_file(args,output_files)

# run main program
if __name__ == '__main__':
    main()
