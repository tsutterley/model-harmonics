#!/usr/bin/env python
u"""
gldas_monthly_harmonics.py
Written by Tyler Sutterley (12/2022)

Reads monthly GLDAS total water storage anomalies and converts to
    spherical harmonic coefficients

Processes as described on the GRACE Tellus Website:
    Data from the Noah 2.7.1 land hydrology model in the Global Land
    Data Assimilation System (GLDAS). The GLDAS system is described
    in the article by Rodell et al (2004). The input data for our
    processing was downloaded from the Goddard Space Flight Center DISC.
The mapped data available at this site is integrated total water content,
    obtained from the GLDAS output by summing the layers:
        Snow Fall water equivalent [kg/m^2]
        Total canopy water storage [kg/m^2]
        Soil Moisture [kg/m^2] (CLM: 10 layers, NOAH: 4 layers)
            CLM:
                0.000 - 0.018 m
                0.018 - 0.045 m
                0.045 - 0.091 m
                0.091 - 0.166 m
                0.166 - 0.289 m
                0.289 - 0.493 m
                0.493 - 0.829 m
                0.829 - 1.383 m
                1.383 - 2.296 m
                2.296 - 3.433 m
            NOAH:
                0.0 - 0.1 m
                0.1 - 0.4 m
                0.4 - 1.0 m
                1.0 - 2.0 m

    Time-averaged grid from a set yearly range subtracted from individual grids.

CALLING SEQUENCE:
    python gldas_monthly_harmonics.py --lmax 60 --format netCDF4 NOAH

INPUTS:
    GLDAS land surface model
        CLM: Common Land Model (CLM)
        CLSM: Catchment Land Surface Model (CLSM)
        MOS: Mosaic model
        NOAH: Noah model
        VIC: Variable Infiltration Capacity (VIC) model

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to run
    -v X, --version X: GLDAS model version
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
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
        can use a custom set of masks to reduce terrestrial water storage
        can use variable loglevels for verbose output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 10/2021: using python logging for handling verbose output
        use output harmonic file wrapper routine to write to file
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 02/2021: include GLDAS MOD44W land mask modified for HYMAP
        replaced numpy bool to prevent deprecation warning
    Updated 01/2021: harmonics object output from gen_stokes.py
    Updated 12/2020: added more love number options
        set spatial variables for both 025 and 10 cases
        using utilities from time module. added maximum harmonic order option
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
        use utilities to define path to load love numbers file
    Updated 06/2020: using spatial data class for input and output operations
    Updated 04/2020: updates to reading load love numbers
        using harmonics class for outputting data to file
    Updated 10/2019: changing Y/N flags to True/False
    Updated 07/2019: output index and date files in separate loop for all files
    Updated 09/2018: use permafrost index from permafrost_gldas_mask.py
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: include Svalbard and Iceland in combined land mask
        output combined land mask to netCDF4 file for validation
    Updated 01/2018: using getopt to set parameters
    Updated 08/2017: convert from geodetic coordinates to geocentric
    Updated 06/2016: updated to use __future__ print function
    Updated 05/2016: complete rewrite of program
    Updated 02/2014: quick code updates
    Written 04/2013
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

# GLDAS models
gldas_products = {}
gldas_products['CLM'] = 'GLDAS Common Land Model (CLM)'
gldas_products['CLSM'] = 'GLDAS Catchment Land Surface Model (CLSM)'
gldas_products['MOS'] = 'GLDAS Mosaic model'
gldas_products['NOAH'] = 'GLDAS Noah model'
gldas_products['VIC'] = 'GLDAS Variable Infiltration Capacity (VIC) model'

# PURPOSE: convert GLDAS terrestrial water storage data to spherical harmonics
def gldas_monthly_harmonics(ddir, MODEL, YEARS,
    SPACING=None,
    VERSION=None,
    MASKS=None,
    LMAX=0,
    MMAX=None,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MODE=0o775):

    # Version flags
    V1,V2 = (f'_V{VERSION}','') if (VERSION == '1') else ('',f'.{VERSION}')
    # use GLDAS monthly products
    TEMPORAL = 'M'
    # subdirectory for model monthly products at spacing for version
    subdir = f'GLDAS_{MODEL}{SPACING}_{TEMPORAL}{V2}'
    # Creating output subdirectory if it doesn't exist
    output_sub = f'GLDAS_{MODEL}{SPACING}{V1}_TWC_CLM_L{LMAX:d}'
    if (not os.access(os.path.join(ddir,output_sub), os.F_OK)):
        os.makedirs(os.path.join(ddir,output_sub),MODE)
    # attributes for output files
    attributes = {}
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # parameters for each grid spacing
    if (SPACING == '025'):
        nlon, nlat = (1440, 600)
        dlon,dlat = (0.25, 0.25)
        extent = [-179.875, 179.875, -59.875, 89.875]
    elif (SPACING == '10'):
        nlon,nlat = (360, 150)
        dlon,dlat = (1.0, 1.0)
        extent = [-179.5, 179.5, -59.5, 89.5]

    # GLDAS MOD44W land mask modified for HYMAP
    landmask_file = f'GLDASp5_landmask_{SPACING}d.nc4'
    with netCDF4.Dataset(os.path.join(ddir,landmask_file),'r') as fileID:
        GLDAS_mask = fileID.variables['GLDAS_mask'][:].squeeze()
        glon = fileID.variables['lon'][:].copy()
        glat = fileID.variables['lat'][:].copy()
    # create mesh grid of latitude and longitude
    gridlon,gridlat = np.meshgrid(glon,glat)
    # create combined mask
    combined_mask = np.logical_not(GLDAS_mask)

    if MASKS:
        # read masks for reducing regions before converting to harmonics
        for mask_file in MASKS:
            fileID = netCDF4.Dataset(mask_file,'r')
            combined_mask |= fileID.variables['mask'][:].astype(bool)
            fileID.close()
    else:
        # use default masks for reducing regions before converting to harmonics
        # mask combining vegetation index, permafrost index and Arctic mask
        # read vegetation index file
        vegetation_file = f'modmodis_domveg20_{SPACING}.nc'
        with netCDF4.Dataset(os.path.join(ddir,vegetation_file),'r') as fileID:
            vegetation_index = fileID.variables['index'][:].copy()
        # 0: missing value
        # 13: Urban and Built-Up
        # 15: Snow and Ice
        # 17: Ocean
        # 18: Wooded Tundra
        # 19: Mixed Tundra
        # 20: Bare Ground Tundra
        for invalid_keys in (0,13,15,17,18,19,20):
            combined_mask |= (vegetation_index == invalid_keys)
        # read Permafrost index file
        permafrost_file = f'permafrost_mod44w_{SPACING}.nc'
        with netCDF4.Dataset(os.path.join(ddir,permafrost_file),'r') as fileID:
            permafrost_index = fileID.variables['mask'][:]
        # 1: Continuous Permafrost
        # 2: Discontinuous Permafrost
        # 3: Isolated Permafrost
        # 4: Sporadic Permafrost
        # 5: Glaciated Area
        for invalid_keys in (1,5):
            combined_mask |= (permafrost_index == invalid_keys)
        # read Arctic mask file
        arctic_file = f'arcticmask_mod44w_{SPACING}.nc'
        with netCDF4.Dataset(os.path.join(ddir,arctic_file),'r') as fileID:
            arctic_mask = fileID.variables['mask'][:].astype(bool)
        # arctic mask
        combined_mask |= arctic_mask[:,:]

    # Earth Parameters
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

    # find input terrestrial water storage files
    regex_years = r'\d+' if (YEARS is None) else r'|'.join(map(str,YEARS))
    args = (MODEL, SPACING, regex_years, suffix[DATAFORM])
    rx = re.compile(r'GLDAS_{0}{1}_TWC_({2})_(\d+)\.{3}$'.format(*args))
    FILES = sorted([fi for fi in os.listdir(os.path.join(ddir,subdir))
        if rx.match(fi)])

    # for each input file
    for t,fi in enumerate(FILES[:-1]):
        # extract year and month from file
        YY,MM = np.array(rx.findall(fi).pop(), dtype=np.float64)

        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            M1 = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,subdir,fi))
            M2 = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,nlon=nlon,
                extent=extent).from_ascii(os.path.join(ddir,subdir,FILES[t+1]))
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            M1 = gravtk.spatial().from_netCDF4(os.path.join(ddir,subdir,fi))
            M2 = gravtk.spatial().from_netCDF4(os.path.join(ddir,subdir,FILES[t+1]))
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            M1 = gravtk.spatial().from_HDF5(os.path.join(ddir,subdir,fi))
            M2 = gravtk.spatial().from_HDF5(os.path.join(ddir,subdir,FILES[t+1]))

        # replace fill value points and certain vegetation types with 0
        M1.replace_invalid(0.0, mask=combined_mask)
        M2.replace_invalid(0.0, mask=combined_mask)
        # calculate 2-month moving average
        # weighting by number of days in each month
        dpm = gravtk.time.calendar_days(int(YY))
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        MASS = (dpm[t % 12]*M1.data + dpm[(t+1) % 12]*M2.data)/W

        # convert to spherical harmonics
        gldas_Ylms = gravtk.gen_stokes(MASS, glon, latitude_geocentric[:,0],
            LMAX=LMAX, MMAX=MMAX, PLM=PLM, LOVE=LOVE)
        # calculate date information
        gldas_Ylms.time, = gravtk.time.convert_calendar_decimal(YY,MM)
        # calculate GRACE/GRACE-FO month
        gldas_Ylms.month = gravtk.time.calendar_to_grace(YY,MM)

        # output spherical harmonic data file
        args=(MODEL,SPACING,LMAX,order_str,gldas_Ylms.month,suffix[DATAFORM])
        FILE='GLDAS_{0}{1}_TWC_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)
        gldas_Ylms.to_file(os.path.join(ddir,output_sub,FILE),
            format=DATAFORM, **attributes)
        # change the permissions mode of the output file to MODE
        os.chmod(os.path.join(ddir,output_sub,FILE),MODE)

    # Output date ascii file
    output_date_file = f'GLDAS_{MODEL}{SPACING}_TWC_DATES.txt'
    fid1 = open(os.path.join(ddir,output_sub,output_date_file),
        mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = 'index.txt'
    fid2 = open(os.path.join(ddir,output_sub,output_index_file),
        mode='w', encoding='utf8')
    # find all available output files
    args = (MODEL, SPACING, LMAX, order_str, suffix[DATAFORM])
    output_regex=r'GLDAS_{0}{1}_TWC_CLM_L{2:d}{3}_([-]?\d+).{4}'.format(*args)
    # find all output ECCO OBP harmonic files (not just ones created in run)
    output_files=[fi for fi in os.listdir(os.path.join(ddir,output_sub))
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
        description="""Reads monthly GLDAS total water storage anomalies
            and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+', choices=gldas_products.keys(),
        help='GLDAS land surface model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    # GLDAS model version
    parser.add_argument('--version','-v',
        type=str, default='2.1',
        help='GLDAS model version')
    # model spatial resolution
    # 10: 1.0 degrees latitude/longitude
    # 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to run')
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

    # for each GLDAS model
    for MODEL in args.model:
        # run program
        gldas_monthly_harmonics(args.directory, MODEL, args.year,
            VERSION=args.version,
            SPACING=args.spacing,
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
