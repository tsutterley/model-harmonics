#!/usr/bin/env python
u"""
merra_smb_harmonics.py
Written by Tyler Sutterley (12/2022)
Reads monthly MERRA-2 surface mass balance anomalies and
    converts to spherical harmonic coefficients

https://disc.gsfc.nasa.gov/information/documents?title=
Records%20of%20MERRA-2%20Data%20Reprocessing%20and%20Service%20Changes

INPUTS:
    SMB: Surface Mass Balance
    ACCUM: Snowfall accumulation
    PRECIP: Total Precipitation
    RAINFALL: Total Rainfall
    SUBLIM: Evaporation and Sublimation
    RUNOFF: Meltwater Runoff

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
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: convert to mid-month averages to correspond with GRACE
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
        added checks for previous versions of reprocessed files
        use output harmonic file wrapper routine to write to file
        add more derived products
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 08/2021: set all points to valid if not using masks
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 02/2021: can use multiple mask files to create a combined solution
        replaced numpy bool to prevent deprecation warning
    Updated 01/2021: added more love number options
        set spatial variables for both 025 and 10 cases
        using utilities from time module. added maximum harmonic order option
        harmonics object output from gen_stokes.py
    Updated 10/2020: use argparse to set command line parameters
    Updated 10/2019: changing Y/N flags to True/False
    Updated 08/2017: convert from geodetic coordinates to geocentric
    Updated 01/2017: can output different data products (SMB, PRECIP, RUNOFF)
    Updated 11/2016: changes to shapefile read for Antarctica
    Written 11/2016
"""
from __future__ import print_function

import sys
import os
import re
import copy
import logging
import netCDF4
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read Merra-2 cumulative data and convert to spherical harmonics
def merra_smb_harmonics(ddir, PRODUCT, YEARS, RANGE=None, REGION=None,
    MASKS=None, LMAX=0, MMAX=None, LOVE_NUMBERS=0, REFERENCE=None,
    DATAFORM=None, MODE=0o775):

    # setup subdirectories
    cumul_sub = f'{PRODUCT}.5.12.4.CUMUL.{RANGE[0]:d}.{RANGE[1]:d}'
    # Creating output subdirectory if it doesn't exist
    prefix = f'{REGION}_' if REGION else ''
    output_sub = f'{prefix}{PRODUCT}_5.12.4_CUMUL_CLM_L{LMAX:d}'
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    # titles for each output data product
    merra_products = {}
    merra_products['SMB'] = 'MERRA-2 Surface Mass Balance'
    merra_products['ACCUM'] = 'MERRA-2 Snowfall accumulation'
    merra_products['PRECIP'] = 'MERRA-2 Precipitation'
    merra_products['RAINFALL'] = 'MERRA-2 Rainfall'
    merra_products['SUBLIM'] = 'MERRA-2 Evaporation and Sublimation'
    merra_products['RUNOFF'] = 'MERRA-2 Meltwater Runoff'
    # source of each output data product
    merra_sources = {}
    merra_sources['SMB'] = ['PRECCU','PRECLS','PRECSN','EVAP','RUNOFF','WESNSC']
    merra_sources['ACCUM'] = ['PRECSN','EVAP']
    merra_sources['PRECIP'] = ['PRECCU','PRECLS','PRECSN']
    merra_sources['RAINFALL'] = ['PRECCU','PRECLS']
    merra_sources['SUBLIM'] = ['EVAP','WESNSC']
    merra_sources['RUNOFF'] = ['RUNOFF']
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # attributes for output files
    attributes = {}
    attributes['title'] = copy.copy(merra_products[PRODUCT])
    attributes['source'] = ', '.join(merra_sources[PRODUCT])
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''

    # output dimensions and extents
    nlat,nlon = (361,576)
    extent = [-180.0,179.375,-90.0,90.0]
    # grid spacing
    dlon,dlat = (0.625,0.5)
    # latitude and longitude
    glon = np.arange(extent[0],extent[1]+dlon,dlon)
    glat = np.arange(extent[2],extent[3]+dlat,dlat)
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
    # semimajor axis of ellipsoid [cm]
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

    # find input files from merra_smb_cumulative.py
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (PRODUCT, regex_years, suffix[DATAFORM])
    regex_pattern = r'MERRA2_(\d+).tavgM_2d_{0}_cumul_Nx.(({1})(\d{{2}})).{2}$'
    rx = re.compile(regex_pattern.format(*args), re.VERBOSE)
    # will be out of order for 2020 due to September reprocessing
    FILES = sorted([fi for fi in os.listdir(os.path.join(ddir,cumul_sub))
        if rx.match(fi)])
    # sort files by month
    indices = np.argsort([rx.match(f1).group(2) for f1 in FILES])
    FILES = [FILES[indice] for indice in indices]
    # remove files that needed to be reprocessed
    INVALID = []
    INVALID.append('MERRA2_400.tavgM_2d_{0}_cumul_Nx.202009.{2}'.format(*args))
    if (set(INVALID) & set(FILES)):
        logging.warning("Reprocessed file found in list")
        FILES = sorted(set(FILES) - set(INVALID))

    # for each input file
    for t,fi in enumerate(FILES[:-1]):
        # extract parameters from input flux file
        MOD,_,YY,MM = rx.findall(fi).pop()
        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            M1 = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,cumul_sub,fi))
            M2 = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,cumul_sub,FILES[t+1]))
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            M1 = gravtk.spatial().from_netCDF4(os.path.join(ddir,cumul_sub,fi),
                varname=PRODUCT)
            M2 = gravtk.spatial().from_netCDF4(os.path.join(ddir,cumul_sub,FILES[t+1]),
                varname=PRODUCT)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            M1 = gravtk.spatial().from_HDF5(os.path.join(ddir,cumul_sub,fi),
                varname=PRODUCT)
            M2 = gravtk.spatial().from_HDF5(os.path.join(ddir,cumul_sub,FILES[t+1]),
                varname=PRODUCT)

        # if reducing to a region of interest before converting to harmonics
        if np.any(input_mask):
            # replace fill value points and masked points with 0
            M1.replace_invalid(0.0, mask=np.logical_not(input_mask))
            M2.replace_invalid(0.0, mask=np.logical_not(input_mask))
        else:
            # replace fill value points points with 0
            M1.replace_invalid(0.0)
            M2.replace_invalid(0.0)

        # calculate 2-month moving average
        # weighting by number of days in each month
        dpm = gravtk.time.calendar_days(int(YY))
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        MASS = (dpm[t % 12]*M1.data + dpm[(t+1) % 12]*M2.data)/W
        # convert to spherical harmonics from mm w.e.
        merra_Ylms = gravtk.gen_stokes(MASS, glon, latitude_geocentric[:,0],
            LMAX=LMAX, MMAX=MMAX, UNITS=3, PLM=PLM, LOVE=LOVE)
        # copy date information
        merra_Ylms.time = np.mean([M1.time, M2.time])
        # calculate GRACE/GRACE-FO month
        merra_Ylms.month = gravtk.time.calendar_to_grace(
            np.float64(YY), np.float64(MM))
        # output spherical harmonic data file
        args=(MOD,PRODUCT,LMAX,order_str,merra_Ylms.month,suffix[DATAFORM])
        FILE='MERRA2_{0}_tavgM_2d_{1}_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)
        merra_Ylms.to_file(os.path.join(ddir,output_sub,FILE),
            format=DATAFORM, **attributes)
        # change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    # Output date ascii file
    output_date_file = f'MERRA2_{PRODUCT}_DATES.txt'
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),
        mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),
        mode='w', encoding='utf8')
    # find all available output files
    args = (PRODUCT, LMAX, order_str, suffix[DATAFORM])
    output_pattern = r'MERRA2_(\d+)_tavgM_2d_{0}_CLM_L{1:d}{2}_([-]?\d+).{3}'
    output_regex = re.compile(output_pattern.format(*args), re.VERBOSE)
    # find all output ECCO OBP harmonic files (not just ones created in run)
    output_files=[fi for fi in os.listdir(os.path.join(ddir,output_sub))
        if re.match(output_regex,fi)]
    for fi in sorted(output_files):
        # extract GRACE month
        MOD,grace_month=np.array(re.findall(output_regex,fi).pop(),dtype=np.int64)
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
        description="""Reads monthly MERRA-2 surface mass balance
            anomalies and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['SMB','ACCUM','PRECIP','RAINFALL','SUBLIM','RUNOFF']
    parser.add_argument('product',
        type=str, nargs='+', choices=choices,
        help='MERRA-2 derived product')
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

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program for each input product
    for PRODUCT in args.product:
        # run program
        merra_smb_harmonics(args.directory, PRODUCT, args.year,
            RANGE=args.mean,
            REGION=args.region,
            MASKS=args.mask,
            LMAX=args.lmax,
            MMAX=args.mmax,
            LOVE_NUMBERS=args.love,
            REFERENCE=args.reference,
            DATAFORM=args.format,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
