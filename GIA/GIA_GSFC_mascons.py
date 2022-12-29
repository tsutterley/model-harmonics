#!/usr/bin/env python
u"""
GIA_GSFC_mascons.py
Written by Tyler Sutterley (12/2022)
Calculates GIA equivalent water height corrections for GSFC mascons at the
    central points of each mascon

INPUTS:
    grace_file: GSFC GRACE/GRACE-FO mascon file to read

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    -G X, --GIA X: GIA model type to read and output
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
    -l X, --lmax X: maximum degree of spherical harmonics
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -V, --verbose: verbose output of processing run
    -M X, --mode X: permissions mode of output GIA mascon file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Pythonic interface to the HDF5 binary data format
        (http://www.h5py.org/)
    netCDF4: Python interface to the netCDF C library
         (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files
    read_GIA_model.py: reads harmonics for a glacial isostatic adjustment model
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    clenshaw_summation.py: calculate spatial field from spherical harmonics
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units

REFERENCES:
    Holmes and Featherstone, "A Unified Approach to the Clenshaw Summation and
        the Recursive Computation of Very High Degree and Order Normalised
        Associated Legendre Functions", Journal of Geodesy (2002)
        http://dx.doi.org/10.1007/s00190-002-0216-2
    Tscherning and Poder, "Some Geodetic Applications of Clenshaw Summation",
        Bollettino di Geodesia e Scienze (1982)

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use full data citation from GIA import function
    Updated 04/2022: use wrapper function for reading load Love numbers
        use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 03/2021: added more love number options
    Updated 12/2020: using argparse to set parameters
    Updated 04/2020: updates to reading load love numbers
    Updated 09/2019: added options for AW13-ICE6G, Caron et al and ICE6G-D
        include lon_span and lat_span in output mascon file
    Written 07/2018
"""
from __future__ import print_function

import sys
import os
import re
import time
import h5py
import logging
import argparse
import traceback
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: keep track of threads
def info(args):
    logging.info(os.path.basename(sys.argv[0]))
    logging.info(args)
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: calculate GIA corrections at locations of GSFC mascons
def GIA_GSFC_mascons(grace_file, GIA=None, GIA_FILE=None, LMAX=0,
    LOVE_NUMBERS=None, REFERENCE=None, MODE=None):

    # set the GRACE directory
    grace_dir = os.path.dirname(grace_file)
    # extract mascon release and version from filename
    rx = re.compile(r'(GSFC\.glb|gsfc.glb_)\.(\d{4})(\d{2})\_'
        r'(\d{4})(\d{2})\_(.*?)(\_.*?)?.h5',re.VERBOSE)
    PROC,SY,SM,EY,EM,VERSION,AUX = rx.findall(grace_file).pop()
    # read the HDF5 file
    output_data = {}
    logging.info(f'{grace_file} -->')
    with h5py.File(grace_file,'r') as fileID:
        for key in ['lat_center','lon_center','lat_span','lon_span']:
            output_data[key] = fileID['mascon'][key][:].flatten()

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # input GIA spherical harmonic datafiles
    GIA_Ylms = gravtk.read_GIA_model(GIA_FILE, GIA=GIA, LMAX=LMAX)
    # Converting GIA rates to cm w.e. for mascon coordinates
    gia_output = gravtk.clenshaw_summation(GIA_Ylms['clm'], GIA_Ylms['slm'],
        output_data['lon_center'], output_data['lat_center'],
        RAD=0, UNITS=1, LMAX=LMAX, LOVE=LOVE)
    output_data['gia'] = gia_output.astype(np.float64)

    # output to file
    FILE = f'GIA_{GIA_Ylms["title"]}_L{LMAX:d}_GSFC_mascons.h5'
    logging.info('\t{0}'.format(os.path.join(grace_dir, FILE)))
    HDF5_GSFC_mascons(output_data, GIA_Ylms, VERSION=VERSION,
        FILENAME=os.path.join(grace_dir, FILE), CLOBBER=True)
    # change the permissions mode
    os.chmod(os.path.join(grace_dir, FILE), MODE)

# PURPOSE: outputting the interpolated data to HDF5
def HDF5_GSFC_mascons(output_data, GIA,
    VERSION='', FILENAME='', CLOBBER=False):

    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)
    # create sub-groups within HDF5 file
    fileID.create_group('mascon')
    fileID.create_group('correction')

    # Dimensions of parameters
    nmas, = output_data['gia'].shape

    # HDF5 file attributes
    attrib = {}
    # latitude
    attrib['lat_center'] = {}
    attrib['lat_center']['long_name'] = 'Latitude_of_center_of_mascon'
    attrib['lat_center']['units'] = 'Degrees_North'
    attrib['lat_span'] = {}
    attrib['lat_span']['long_name'] = 'Mascon_latitude_central_angle'
    attrib['lat_span']['units'] = 'Degrees'
    # longitude
    attrib['lon_center'] = {}
    attrib['lon_center']['long_name'] = 'Longitude_of_center_of_mascon'
    attrib['lon_center']['units'] = 'Degrees_East'
    attrib['lon_span'] = {}
    attrib['lon_span']['long_name'] = 'Mascon_longitude_central_angle'
    attrib['lon_span']['units'] = 'Degrees'
    # GIA correction
    attrib['gia'] = {}
    attrib['gia']['long_name'] = 'Equivalent_water_height'
    attrib['gia']['model'] = GIA['title']
    attrib['gia']['description'] = ('Equivalent water thickness rate due '
        'to Glacial Isostatic Adjustment (GIA)')
    attrib['gia']['source'] = GIA['citation']
    attrib['gia']['reference'] = GIA['reference']
    attrib['gia']['coordinates'] = '/mascon/lat_center /mascon/lon_center'
    attrib['gia']['units'] = 'cm/yr'
    # groups for each key
    groups = dict(lat_center='mascon', lon_center='mascon', lat_span='mascon',
        lon_span='mascon', gia='correction')

    # create HDF5 records
    h5 = {}
    for key,val in output_data.items():
        # Defining the HDF5 dataset variables
        h5[key] = fileID.create_dataset(f'{groups[key]}/{key}',
            (nmas,), data=val, dtype=val.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in attrib[key].items():
            h5[key].attrs[att_name] = att_val
        # attach dimensions
        h5[key].dims[0].label = 'RECORD_SIZE'

    # HDF5 file title
    fileID.attrs['featureType'] = 'timeSeries'
    fileID.attrs['title'] = 'GIA_Correction'
    fileID.attrs['summary'] = ('Glacial_Isostatic_Adjustment_Corrections_'
        'for_NASA_Goddard_Space_Flight_Center_(GSFC)_mascons.')
    fileID.attrs['date_created'] = time.strftime('%Y-%m-%d',time.localtime())
    fileID.attrs['version'] = VERSION
    # add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = output_data['lat_center'].min()
    fileID.attrs['geospatial_lat_max'] = output_data['lat_center'].max()
    fileID.attrs['geospatial_lon_min'] = output_data['lon_center'].min()
    fileID.attrs['geospatial_lon_max'] = output_data['lon_center'].max()
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    # add software information
    fileID.attrs['software_reference'] = mdlhmc.version.project_name
    fileID.attrs['software_version'] = mdlhmc.version.full_version
    # Closing the HDF5 file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads dates of GSFC GRACE mascon data files and assigns
            the month number
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='GSFC GRACE mascon file to run')
    # GIA model type list
    models = {}
    models['IJ05-R2'] = 'Ivins R2 GIA Models'
    models['W12a'] = 'Whitehouse GIA Models'
    models['SM09'] = 'Simpson/Milne GIA Models'
    models['ICE6G'] = 'ICE-6G GIA Models'
    models['Wu10'] = 'Wu (2010) GIA Correction'
    models['AW13-ICE6G'] = 'Geruo A ICE-6G GIA Models'
    models['AW13-IJ05'] = 'Geruo A IJ05-R2 GIA Models'
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
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GIA file to read')
    # degree of truncation
    parser.add_argument('--lmax','-l',
        type=int, default=120,
        help='Maximum degree of spherical harmonics')
    # different treatments of the load Love numbers
    # 0: Han and Wahr (1995) values from PREM
    # 1: Gegout (2005) values from PREM
    # 2: Wang et al. (2012) values from PREM
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2],
        help='Treatment of the Load Love numbers')
    # option for setting reference frame for gravitational load love number
    # reference frame options (CF, CM, CE)
    parser.add_argument('--reference',
        type=str.upper, default='CF', choices=['CF','CM','CE'],
        help='Reference frame for load Love numbers')
    # print information about processing run
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # calculate GIA rates for each GSFC GRACE mascon data
    for grace_file in args.infile:
        # try to run the analysis with listed parameters
        try:
            info(args)
            GIA_GSFC_mascons(grace_file,
                GIA=args.gia,
                GIA_FILE=args.gia_file,
                LMAX=args.lmax,
                LOVE_NUMBERS=args.love,
                REFERENCE=args.reference,
                MODE=args.mode)
        except Exception as e:
            # if there has been an error exception
            # print the type, value, and stack trace of the
            # current exception being handled
            logging.critical(f'process id {os.getpid():d} failed')
            logging.error(traceback.format_exc())

# run main program
if __name__ == '__main__':
    main()
