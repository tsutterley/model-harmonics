#!/usr/bin/env python
u"""
harmonic_operators.py
Written by Tyler Sutterley (02/2021)
Performs basic operations on spherical harmonic files

CALLING SEQUENCE:
    python harmonic_operators.py --operation add infile1 infile2 outfile

INPUTS:
    path to input harmonic files
    path to output harmonic file

COMMAND LINE OPTIONS:
    -O X, --operation X: Operation to run
        add
        subtract
        multiply
        divide
        mean
        destripe
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
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
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
        destripe_harmonics.py: calculates the decorrelation (destriping) filter
            and filters the GRACE/GRACE-FO coefficients for striping errors
        ncdf_read_stokes.py: reads spherical harmonic netcdf files
        ncdf_stokes.py: writes output spherical harmonic data to netcdf
        hdf5_read_stokes.py: reads spherical harmonic HDF5 files
        hdf5_stokes.py: writes output spherical harmonic data to HDF5

UPDATE HISTORY:
    Updated 02/2021: added options to truncate output to a degree or order
        add options to read from individual index files
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import argparse
import numpy as np
from gravity_toolkit.harmonics import harmonics

#-- PURPOSE: Performs operations on harmonic files
def harmonic_operators(INPUT_FILES, OUTPUT_FILE, OPERATION=None, LMAX=None,
    MMAX=None, DATAFORM=None, DATE=False, VERBOSE=False, MODE=None):

    #-- number of input harmonic files
    n_files = len(INPUT_FILES)
    #-- extend list if a single format was entered for all files
    if len(DATAFORM) < (n_files+1):
        DATAFORM = DATAFORM*(n_files+1)
    #-- verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY,MODE,exist_ok=True)

    #-- read each input file
    dinput = [None]*n_files
    for i,fi in enumerate(INPUT_FILES):
        #-- read spherical harmonics file in data format
        if DATAFORM[i] in ('ascii','netCDF4','HDF5'):
            #-- ascii (.txt)
            #-- netCDF4 (.nc)
            #-- HDF5 (.H5)
            dinput[i] = harmonics().from_file(fi,format=DATAFORM[i],
                date=DATE, verbose=VERBOSE)
        elif DATAFORM[i] in ('index-ascii','index-netCDF4','index-HDF5'):
            #-- read from index file
            _,dataform = DATAFORM[i].split('-')
            dinput[i] = harmonics().from_index(fi,format=dataform,date=DATE)

    #-- operate on input files
    if (OPERATION == 'add'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            output = output.add(dinput[i])
    elif (OPERATION == 'subtract'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            #-- perform operation
            output = output.subtract(dinput[i+1])
    elif (OPERATION == 'multiply'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            #-- perform operation
            output = output.multiply(dinput[i+1])
    elif (OPERATION == 'divide'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            #-- perform operation
            output = output.divide(dinput[i+1])
    elif (OPERATION == 'mean'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            #-- perform operation
            output = output.add(dinput[i])
        #-- convert from total to mean
        output = output.scale(1.0/n_files)
    elif (OPERATION == 'destripe'):
        #-- destripe spherical harmonics
        output = dinput[0].destripe()
    #-- truncate to specified degree and order
    if (LMAX is not None) | (MMAX is not None):
        output.truncate(LMAX, mmax=MMAX)
    #-- copy date variables if specified
    if DATE:
        output.time = np.copy(dinput[0].time)
        output.month = np.copy(dinput[0].month)

    #-- write spherical harmonic file in data format
    if (DATAFORM[-1] == 'ascii'):
        #-- ascii (.txt)
        print(OUTPUT_FILE) if VERBOSE else None
        output.to_ascii(OUTPUT_FILE,date=DATE)
    elif (DATAFORM[-1] == 'netCDF4'):
        #-- netcdf (.nc)
        output.to_netCDF4(OUTPUT_FILE,date=DATE,VERBOSE=VERBOSE,
            TITLE='Output from {0}'.format(os.path.basename(sys.argv[0])))
    elif (DATAFORM[-1] == 'HDF5'):
        #-- HDF5 (.H5)
        output.to_HDF5(OUTPUT_FILE,date=DATE,VERBOSE=VERBOSE,
            TITLE='Output from {0}'.format(os.path.basename(sys.argv[0])))
    #-- change the permissions mode of the output file
    os.chmod(OUTPUT_FILE, MODE)

#-- Main program that calls harmonic_operators()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Performs basic operations on spherical harmonic files
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
        choices=['add','subtract','multiply','divide','mean','destripe'],
        help='Operation to run')
    #-- maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=None,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
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
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- run program
    harmonic_operators(args.infiles, args.outfile[0], OPERATION=args.operation,
        LMAX=args.lmax, MMAX=args.mmax, DATAFORM=args.format, DATE=args.date,
        VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()