#!/usr/bin/env python
u"""
spatial_operators.py
Written by Tyler Sutterley (05/2022)
Performs basic operations on spatial files

CALLING SEQUENCE:
    python spatial_operators.py --operation add infile1 infile2 outfile

INPUTS:
    path to input spatial files
    path to output spatial file

COMMAND LINE OPTIONS:
    -O X, --operation X: Operation to run
        add
        subtract
        multiply
        divide
        mean
        error
        RMS
    -S X, --spacing X: spatial resolution of input data (dlon,dlat)
    -I X, --interval X: input grid interval
        1: (0:360, 90:-90)
        2: (degree spacing/2)
    --header X: number of header rows to skip in input ascii files
    -F X, --format X: Input and output data format
        ascii
        netcdf
        HDF5
    -D, --date: input and output files have date information
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: using python logging for handling verbose output
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: added variance off mean as estimated error
        add options to read from individual index files
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import logging
import argparse
import numpy as np
from gravity_toolkit.spatial import spatial

#-- PURPOSE: Performs operations on spatial files
def spatial_operators(INPUT_FILES, OUTPUT_FILE, OPERATION=None, DDEG=None,
    INTERVAL=None, HEADER=None, DATAFORM=None, DATE=False, MODE=None):

    #-- number of input spatial files
    n_files = len(INPUT_FILES)
    #-- extend list if a single format was entered for all files
    if len(DATAFORM) < (n_files+1):
        DATAFORM = DATAFORM*(n_files+1)
    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- Grid spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    #-- Grid dimensions
    if (INTERVAL == 1):#-- (0:360, 90:-90)
        nlon = np.int64((360.0/dlon)+1.0)
        nlat = np.int64((180.0/dlat)+1.0)
    elif (INTERVAL == 2):#-- degree spacing/2
        nlon = np.int64((360.0/dlon))
        nlat = np.int64((180.0/dlat))

    #-- read each input file
    dinput = [None]*n_files
    for i,fi in enumerate(INPUT_FILES):
        #-- read spatial file in data format
        if DATAFORM[i] in ('ascii','netCDF4','HDF5'):
            #-- ascii (.txt)
            #-- netCDF4 (.nc)
            #-- HDF5 (.H5)
            dinput[i] = spatial(spacing=[dlon,dlat],nlat=nlat,
                nlon=nlon).from_file(fi,format=DATAFORM[i],date=DATE)
        elif DATAFORM[i] in ('index-ascii','index-netCDF4','index-HDF5'):
            #-- read from index file
            _,dataform = DATAFORM[i].split('-')
            dinput[i] = spatial(spacing=[dlon,dlat],nlat=nlat,
                nlon=nlon).from_index(fi,format=dataform,date=DATE)

    #-- operate on input files
    if (OPERATION == 'add'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            output = output.offset(dinput[i].data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i].mask)
    elif (OPERATION == 'subtract'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            #-- perform operation
            output = output.offset(dinput[i+1].scale(-1).data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i+1].mask)
    elif (OPERATION == 'multiply'):
        output = dinput[0].zeros_like().offset(1.0)
        for i in range(n_files):
            #-- perform operation
            output = output.scale(dinput[i].data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i].mask)
    elif (OPERATION == 'divide'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            #-- perform operation
            output = output.scale(dinput[i+1].power(-1).data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i+1].mask)
    elif (OPERATION == 'mean'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            output = output.offset(dinput[i].data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i].mask)
        #-- convert from total to mean
        output = output.scale(1.0/n_files)
    elif (OPERATION == 'error'):
        mean = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            mean = mean.offset(dinput[i].data)
            #-- update mask with values from file
            mean.replace_invalid(mean.fill_value,mask=dinput[i].mask)
        #-- convert from total to mean
        mean = mean.scale(1.0/n_files)
        #-- use variance off mean as estimated error
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            temp = dinput[i].offset(-mean.data)
            output = output.offset(temp.power(2.0).data)
        #-- update mask with values from mean
        output.replace_invalid(output.fill_value,mask=mean.mask)
        #-- calculate RMS of mean differences
        output = output.scale(1.0/(n_files-1.0)).power(0.5)
    elif (OPERATION == 'RMS'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            output = output.offset(dinput[i].power(2.0).data)
            #-- update mask with values from file
            output.replace_invalid(output.fill_value,mask=dinput[i].mask)
        #-- convert from total in quadrature to RMS
        output = output.scale(1.0/n_files).power(0.5)

    #-- write spatial file in data format
    if (DATAFORM[-1] == 'ascii'):
        #-- ascii (.txt)
        output.to_ascii(OUTPUT_FILE,date=DATE)
    elif (DATAFORM[-1] == 'netCDF4'):
        #-- netcdf (.nc)
        attr = dinput[0].attributes['data']
        output.to_netCDF4(OUTPUT_FILE,date=DATE,
            units=attr['units'],longname=attr['long_name'],
            title='Output from {0}'.format(os.path.basename(sys.argv[0])))
    elif (DATAFORM[-1] == 'HDF5'):
        #-- HDF5 (.H5)
        attr = dinput[0].attributes['data']
        output.to_HDF5(OUTPUT_FILE,date=DATE,
            units=attr['units'],longname=attr['long_name'],
            title='Output from {0}'.format(os.path.basename(sys.argv[0])))
    #-- change the permissions mode of the output file
    os.chmod(OUTPUT_FILE, MODE)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Performs basic operations on spatial files
            """
    )
    #-- command line options
    #-- input and output file
    parser.add_argument('infiles',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Input files')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs=1,
        help='Output file')
    #-- operation to run
    parser.add_argument('--operation','-O',
        metavar='OPERATION', type=str, required=True,
        choices=['add','subtract','multiply','divide','mean','error','RMS'],
        help='Operation to run')
    #-- output grid parameters
    parser.add_argument('--spacing','-S',
        type=float, nargs='+', default=[0.5,0.5], metavar=('dlon','dlat'),
        help='Spatial resolution of output data')
    parser.add_argument('--interval','-I',
        type=int, default=2, choices=[1,2,3],
        help='Output grid interval (1: global, 2: centered global)')
    #-- ascii parameters
    parser.add_argument('--header',
        type=int,
        help='Number of header rows to skip in input ascii files')
    #-- input and output data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    parser.add_argument('--format','-F',
        metavar='FORMAT', type=str, nargs='+',
        default=['netCDF4'], choices=choices,
        help='Input and output data format')
    #-- Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    #-- run program
    spatial_operators(args.infiles, args.outfile[0], OPERATION=args.operation,
        DDEG=args.spacing, INTERVAL=args.interval, HEADER=args.header,
        DATAFORM=args.format, DATE=args.date, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()