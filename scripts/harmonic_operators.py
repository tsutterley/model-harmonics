#!/usr/bin/env python
u"""
harmonic_operators.py
Written by Tyler Sutterley (12/2022)
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
        error
        RMS
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

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
        updated GIA reader to be an inheritance of harmonics
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: can read from GIA models for merging or correcting
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: using python logging for handling verbose output
    Updated 08/2021: added variance off mean as estimated error
    Updated 02/2021: added options to truncate output to a degree or order
        add options to read from individual index files
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import logging
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: Performs operations on harmonic files
def harmonic_operators(INPUT_FILES, OUTPUT_FILE, OPERATION=None, LMAX=None,
    MMAX=None, DATAFORM=None, DATE=False, MODE=None):

    # number of input harmonic files
    n_files = len(INPUT_FILES)
    # extend list if a single format was entered for all files
    if len(DATAFORM) < (n_files+1):
        DATAFORM = DATAFORM*(n_files+1)
    # verify that output directory exists
    DIRECTORY = os.path.abspath(os.path.dirname(OUTPUT_FILE))
    if not os.access(DIRECTORY, os.F_OK):
        os.makedirs(DIRECTORY, mode=MODE, exist_ok=True)

    # list of available GIA Models
    GIA = ['IJ05-R2', 'W12a', 'SM09', 'Wu10',
        'AW13-ICE6G', 'AW13-IJ05', 'Caron', 'ICE6G-D']
    # read each input file
    dinput = [None]*n_files
    for i,fi in enumerate(INPUT_FILES):
        # read spherical harmonics file in data format
        if DATAFORM[i] in ('ascii', 'netCDF4', 'HDF5'):
            # ascii (.txt)
            # netCDF4 (.nc)
            # HDF5 (.H5)
            dinput[i] = gravtk.harmonics().from_file(fi,
                format=DATAFORM[i], date=DATE)
        elif DATAFORM[i] in ('index-ascii', 'index-netCDF4', 'index-HDF5'):
            # read from index file
            _,dataform = DATAFORM[i].split('-')
            dinput[i] = gravtk.harmonics().from_index(fi,
                format=dataform, date=DATE)
        elif (DATAFORM[i] in GIA) and DATE:
            # read from GIA file and calculate drift
            temp = gravtk.gia().from_GIA(fi, GIA=DATAFORM[i])
            dinput[i] = temp.drift(dinput[0].time)
        elif DATAFORM[i] in GIA:
            dinput[i] = gravtk.gia().from_GIA(fi, GIA=DATAFORM[i])

    # operate on input files
    if (OPERATION == 'add'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            # perform operation
            output = output.add(dinput[i])
    elif (OPERATION == 'subtract'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            # perform operation
            output = output.subtract(dinput[i+1])
    elif (OPERATION == 'multiply'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            # perform operation
            output = output.multiply(dinput[i+1])
    elif (OPERATION == 'divide'):
        output = dinput[0].copy()
        for i in range(n_files-1):
            # perform operation
            output = output.divide(dinput[i+1])
    elif (OPERATION == 'mean'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            # perform operation
            output = output.add(dinput[i])
        # convert from total to mean
        output = output.scale(1.0/n_files)
    elif (OPERATION == 'destripe'):
        # destripe spherical harmonics
        output = dinput[0].destripe()
    elif (OPERATION == 'error'):
        mean = dinput[0].zeros_like()
        for i in range(n_files):
            # perform operation
            mean = mean.add(dinput[i])
        # convert from total to mean
        mean = mean.scale(1.0/n_files)
        # use variance off mean as estimated error
        output = dinput[0].zeros_like()
        for i in range(n_files):
            # perform operation
            temp = dinput[i].subtract(mean)
            output = output.add(temp.power(2.0))
        # calculate RMS of mean differences
        output = output.scale(1.0/(n_files-1.0)).power(0.5)
    elif (OPERATION == 'RMS'):
        output = dinput[0].zeros_like()
        for i in range(n_files):
            # perform operation
            output = output.add(dinput[i].power(2.0))
        # convert from total in quadrature to RMS
        output = output.scale(1.0/n_files).power(0.5)
    # truncate to specified degree and order
    if (LMAX is not None) | (MMAX is not None):
        output.truncate(LMAX, mmax=MMAX)
    # copy date variables if specified
    if DATE:
        output.time = np.copy(dinput[0].time)
        output.month = np.copy(dinput[0].month)

    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # write spherical harmonic file in data format
    output.to_file(OUTPUT_FILE, format=DATAFORM[-1],
        date=DATE, **attributes)
    # change the permissions mode of the output file
    os.chmod(OUTPUT_FILE, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Performs basic operations on spherical harmonic files
            """
    )
    # command line options
    # input and output file
    parser.add_argument('infiles',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Input files')
    parser.add_argument('outfile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs=1,
        help='Output file')
    # operation to run
    choices = ['add','subtract','multiply','divide','mean',
        'destripe','error','RMS']
    parser.add_argument('--operation','-O',
        metavar='OPERATION', type=str,
        required=True, choices=choices,
        help='Operation to run')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=None,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
    # input and output data format (ascii, netCDF4, HDF5)
    choices = []
    choices.extend(['ascii','netCDF4','HDF5'])
    choices.extend(['index-ascii','index-netCDF4','index-HDF5'])
    choices.extend(['IJ05-R2','W12a','SM09','Wu10','AW13-ICE6G',
        'AW13-IJ05','Caron','ICE6G-D'])
    parser.add_argument('--format','-F',
        metavar='FORMAT', type=str, nargs='+',
        default=['netCDF4'], choices=choices,
        help='Input and output data format')
    # Input and output files have date information
    parser.add_argument('--date','-D',
        default=False, action='store_true',
        help='Input and output files have date information')
    # print information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    harmonic_operators(args.infiles, args.outfile[0], OPERATION=args.operation,
        LMAX=args.lmax, MMAX=args.mmax, DATAFORM=args.format, DATE=args.date,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()