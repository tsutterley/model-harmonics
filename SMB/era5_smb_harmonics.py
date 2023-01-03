#!/usr/bin/env python
u"""
era5_smb_harmonics.py
Written by Tyler Sutterley (12/2022)
Reads monthly ERA5 surface mass balance anomalies and
    converts to spherical harmonic coefficients

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -m X, --mean X: Year range for mean
    -Y X, --year X: Years to run
    -R X, --region X: region name for subdirectory
    --mask X: netCDF4 mask files for reducing to regions
    -l X, --lmax X: maximum spherical harmonic degree
    -m X, --mmax X: maximum spherical harmonic order
    -n X, --love X: Load Love numbers dataset
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
    --reference X: Reference frame for load love numbers
        CF: Center of Surface Figure (default)
        CM: Center of Mass of Earth System
        CE: Center of Mass of Solid Earth
    -F X, --format X: Input and output data format
        ascii
        netcdf
        HDF5
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
    associated_legendre.py: computes fully-normalized associated Legendre polynomials
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    constants.py: calculate reference parameters for common ellipsoids
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: convert to mid-month averages to correspond with GRACE
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Written 10/2021
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read ERA5 cumulative data and convert to spherical harmonics
def era5_smb_harmonics(ddir, YEARS, RANGE=None, REGION=None,
    MASKS=None, LMAX=0, MMAX=None, LOVE_NUMBERS=0, REFERENCE=None,
    DATAFORM=None, MODE=0o775):

    # setup subdirectories
    cumul_sub = 'ERA5-Cumul-P-E-{0:4d}-{1:4d}'.format(*RANGE)
    # Creating output subdirectory if it doesn't exist
    prefix = f'{REGION}_' if REGION else ''
    output_sub = f'{prefix}ERA5_CUMUL_P-E_CLM_L{LMAX:d}'
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    # output data file format and title
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # attributes for output files
    attributes = {}
    attributes['title'] = 'ERA5 Precipitation minus Evaporation'
    attributes['source'] = ', '.join(['tp','e'])
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''

    # output dimensions and extents
    nlat,nlon = (721,1440)
    extent = [0.0,359.75,-90.0,90.0]
    # grid spacing
    dlon,dlat = (0.25,0.25)
    # latitude and longitude
    glon = np.arange(extent[0],extent[1]+dlon,dlon)
    glat = np.arange(extent[3],extent[2]-dlat,-dlat)
    # create mesh grid of latitude and longitude
    gridlon,gridlat = np.meshgrid(glon,glat)

    # create mask object for reducing data
    if bool(MASKS):
        input_mask = np.zeros((nlat,nlon),dtype=bool)
    else:
        input_mask = np.ones((nlat,nlon),dtype=bool)
    # read masks for reducing regions before converting to harmonics
    for mask_file in MASKS:
        logging.info(mask_file)
        fileID = netCDF4.Dataset(mask_file,'r')
        input_mask |= fileID.variables['mask'][:].astype(bool)
        fileID.close()

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.constants(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude in radians
    latitude_geodetic_rad = np.pi*gridlat/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = N * np.cos(latitude_geodetic_rad) * np.cos(np.pi*gridlon/180.0)
    Y = N * np.cos(latitude_geodetic_rad) * np.sin(np.pi*gridlon/180.0)
    Z = (N * (1.0 - ecc1**2.0)) * np.sin(latitude_geodetic_rad)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi
    # colatitude in radians
    theta = (90.0 - latitude_geocentric[:,0])*np.pi/180.0

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # calculate Legendre polynomials
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # find input files from era5_smb_cumulative.py
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    rx = re.compile(r'ERA5\-Cumul\-P-E\-({0})\.nc$'.format(regex_years))
    FILES = [f for f in os.listdir(os.path.join(ddir,cumul_sub)) if rx.match(f)]

    # create list of yearly ERA5 files
    spatial_list = []
    # for each input file
    for t,fi in enumerate(FILES):
        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            era5_data = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,cumul_sub,fi))
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            era5_data = gravtk.spatial().from_netCDF4(
                os.path.join(ddir,cumul_sub,fi), varname='SMB')
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            era5_data = gravtk.spatial().from_HDF5(
                os.path.join(ddir,cumul_sub,fi), varname='SMB')
        # if reducing to a region of interest before converting to harmonics
        if np.any(input_mask):
            # replace fill value points and masked points with 0
            era5_data.replace_invalid(0.0, mask=np.logical_not(input_mask))
        else:
            # replace fill value points points with 0
            era5_data.replace_invalid(0.0)
        # append to spatial list
        spatial_list.append(era5_data)
    # convert to combined data cube and clear memory from spatial list
    era5_data = gravtk.spatial().from_list(spatial_list, clear=True)
    nlat,nlon,nt = era5_data.shape

    # for each month of data
    for i in range(nt-1):
        # convert data to m w.e.
        M1 = era5_data.index(i).scale(1000.0)
        M2 = era5_data.index(i+1).scale(1000.0)
        # calculate 2-month moving average
        # weighting by number of days in each month
        dpm = gravtk.time.calendar_days(np.floor(M1.time))
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        MASS = (dpm[t % 12]*M1.data + dpm[(t+1) % 12]*M2.data)/W
        # convert to spherical harmonics
        era5_Ylms = gravtk.gen_stokes(MASS, glon, latitude_geocentric[:,0],
            LMAX=LMAX, MMAX=MMAX, UNITS=3, PLM=PLM, LOVE=LOVE)
        # copy date information
        era5_Ylms.time = np.mean([M1.time, M2.time])
        era5_Ylms.month = gravtk.time.calendar_to_grace(
            era5_Ylms.time)
        # output spherical harmonic data file
        args = (LMAX, order_str, era5_Ylms.month, suffix[DATAFORM])
        FILE = 'ERA5_CUMUL_P-E_CLM_L{0:d}{1}_{2:03d}.{3}'.format(*args)
        era5_Ylms.to_file(os.path.join(ddir,output_sub,FILE),
            format=DATAFORM, **attributes)
        # change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    # Output date ascii file
    output_date_file = 'ERA5_SMB_DATES.txt'
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),
        mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),
        mode='w', encoding='utf8')
    # find all available output files
    args = (LMAX, order_str, suffix[DATAFORM])
    output_pattern = r'ERA5_CUMUL_P-E_CLM_L{0:d}{1}_([-]?\d+).{2}'
    output_regex = re.compile(output_pattern.format(*args), re.VERBOSE)
    # find all output ECCO OBP harmonic files (not just ones created in run)
    output_files = [fi for fi in os.listdir(os.path.join(ddir,output_sub))
        if re.match(output_regex,fi)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi),dtype=np.int64)
        YY,MM = gravtk.time.grace_to_calendar(grace_month)
        tdec, = gravtk.time.convert_calendar_decimal(YY, MM)
        # full path to output file
        full_output_file = os.path.join(ddir,output_sub,fi)
        # print date, GRACE month and calendar month to date file
        fid1.write('{0:11.6f} {1:03d} {2:02.0f}\n'.format(tdec,grace_month,MM))
        # print output file to index
        print(full_output_file.replace(os.path.expanduser('~'),'~'), file=fid2)
    # close the date and index files
    fid1.close()
    fid2.close()
    # set the permissions level of the output date and index files to MODE
    os.chmod(os.path.join(ddir,output_sub,output_date_file), MODE)
    os.chmod(os.path.join(ddir,output_sub,output_index_file), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ERA5 surface mass balance
            anomalies and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    # region name for subdirectory
    parser.add_argument('--region','-R',
        type=str, default=None,
        help='Region name for subdirectory')
    # mask file for reducing to regions
    parser.add_argument('--mask',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        nargs='+', default=[],
        help='netCDF4 masks file for reducing to regions')
    # maximum spherical harmonic degree and order
    parser.add_argument('--lmax','-l',
        type=int, default=60,
        help='Maximum spherical harmonic degree')
    parser.add_argument('--mmax','-m',
        type=int, default=None,
        help='Maximum spherical harmonic order')
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
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # print information about each input and output file
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

    # run program with parameters
    era5_smb_harmonics(args.directory, args.year, RANGE=args.mean,
        REGION=args.region, MASKS=args.mask, LMAX=args.lmax, MMAX=args.mmax,
        LOVE_NUMBERS=args.love, REFERENCE=args.reference,
        DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
