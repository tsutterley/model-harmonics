#!/usr/bin/env python
u"""
spatial_operators.py
Written by Tyler Sutterley (02/2021)
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
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5

UPDATE HISTORY:
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import argparse
import numpy as np
from gravity_toolkit.spatial import spatial

#-- PURPOSE: Performs operations on spatial files
def spatial_operators(INPUT_FILES, OUTPUT_FILE, OPERATION=None, DDEG=None,
    INTERVAL=None, HEADER=None, DATAFORM=None, DATE=False, VERBOSE=False,
    MODE=None):

    #-- number of input spatial files
    n_files = len(INPUT_FILES)
    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- Grid spacing
    dlon,dlat = (DDEG,DDEG) if (np.ndim(DDEG) == 0) else (DDEG[0],DDEG[1])
    #-- Grid dimensions
    if (INTERVAL == 1):#-- (0:360, 90:-90)
        nlon = np.int((360.0/dlon)+1.0)
        nlat = np.int((180.0/dlat)+1.0)
    elif (INTERVAL == 2):#-- degree spacing/2
        nlon = np.int((360.0/dlon))
        nlat = np.int((180.0/dlat))

    #-- read each input file
    dinput = [None]*n_files
    for i,fi in enumerate(INPUT_FILES):
        #-- read spatial file in data format
        if (DATAFORM == 'ascii'):
            #-- ascii (.txt)
            dinput[i] = spatial(spacing=[dlon,dlat],nlat=nlat,
                nlon=nlon).from_ascii(fi,header=HEADER,date=DATE)
        elif (DATAFORM == 'netCDF4'):
            #-- netcdf (.nc)
            dinput[i] = spatial().from_netCDF4(fi,date=DATE)
        elif (DATAFORM == 'HDF5'):
            #-- HDF5 (.H5)
            dinput[i] = spatial().from_HDF5(fi,date=DATE)

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
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        output.to_ascii(OUTPUT_FILE,date=DATE,verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        attr = dinput[0].attributes['data']
        output.to_netCDF4(OUTPUT_FILE,date=DATE,verbose=VERBOSE,
            units=attr['units'],longname=attr['long_name'],
            title='Output from {0}'.format(os.path.basename(sys.argv[0])))
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        attr = dinput[0].attributes['data']
        output.to_HDF5(OUTPUT_FILE,date=DATE,verbose=VERBOSE,
            units=attr['units'],longname=attr['long_name'],
            title='Output from {0}'.format(os.path.basename(sys.argv[0])))
    #-- change the permissions mode of the output file
    os.chmod(OUTPUT_FILE, MODE)

#-- Main program that calls spatial_operators()
def main():
    #-- Read the system arguments listed after the program
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
        metavar='OPERATION', type=str,
        choices=['add','subtract','multiply','divide','mean','RMS'],
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
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- run program
    spatial_operators(args.infiles, args.outfile[0], OPERATION=args.operation,
        DDEG=args.spacing, INTERVAL=args.interval, HEADER=args.header,
        DATAFORM=args.format, DATE=args.date, VERBOSE=args.verbose,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()