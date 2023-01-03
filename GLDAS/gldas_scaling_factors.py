#!/usr/bin/env python
u"""
gldas_scaling_factors.py
Written by Tyler Sutterley (12/2022)

Reads monthly GLDAS total water storage anomalies and monthly
    spherical harmonic coefficients
Computes point scaling factors following Landerer and Swenson (2012)
    https://doi.org/10.1029/2011WR011453

CALLING SEQUENCE:
    python gldas_scaling_factors.py --lmax 60 --radius 300 --destripe \
        --format netCDF4 NOAH

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
    -s X, --start X: Starting GRACE/GRACE-FO month
    -e X, --end X: Ending GRACE/GRACE-FO month
    -o X, --missing X: Missing GRACE/GRACE-FO months
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -v X, --version X: GLDAS model version
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
    -R X, --radius X: Gaussian smoothing radius (km)
    -d, --destripe: use a decorrelation filter (destriping filter)
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
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    constants.py: calculate reference parameters for common ellipsoids
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    gauss_weights.py: Computes the Gaussian weights as a function of degree
    harmonic_summation.py: calculates a spatial field from spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    units.py: class for converting GRACE/GRACE-FO Level-2 data to specific units
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 08/2022: create verbose logger within main definition
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 02/2021: output spatial power of original data
        use GRACE/GRACE-FO months to select range of GLDAS data
        include GLDAS MOD44W land mask modified for HYMAP
        replaced numpy bool to prevent deprecation warning
    Updated 01/2021 for public release
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
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import argparse
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

# PURPOSE: read GLDAS terrestrial water storage data and spherical harmonics
# calculate the point scaling factors following Landerer and Swenson (2012)
def gldas_scaling_factors(ddir, MODEL, START_MON, END_MON, MISSING,
    SPACING=None,
    VERSION=None,
    MASKS=None,
    LMAX=0,
    MMAX=None,
    RAD=0,
    DESTRIPE=False,
    LOVE_NUMBERS=0,
    REFERENCE=None,
    DATAFORM=None,
    MODE=0o775):

    # Version flags
    V1,V2 = (f'_V{VERSION}','') if (VERSION == '1') else ('',f'.{VERSION}')
    # use GLDAS monthly products
    TEMPORAL = 'M'
    # subdirectory for model monthly products at spacing for version
    sub1 = f'GLDAS_{MODEL}{SPACING}_{TEMPORAL}{V2}'
    # Creating output subdirectory if it doesn't exist
    sub2 = f'GLDAS_{MODEL}{SPACING}{V1}_TWC_CLM_L{LMAX:d}'
    if (not os.access(os.path.join(ddir,sub2), os.F_OK)):
        os.makedirs(os.path.join(ddir,sub2),MODE)
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # Calculating the Gaussian smoothing for radius RAD
    gw_str = f'_r{RAD:0.0f}km' if (RAD != 0) else ''
    # destriped GRACE/GRACE-FO coefficients
    ds_str = '_FL' if DESTRIPE else ''
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # parameters for each grid spacing
    if (SPACING == '025'):
        nlon,nlat = (1440,600)
        dlon,dlat = (0.25,0.25)
        extent = [-179.875,179.875,-59.875,89.875]
    elif (SPACING == '10'):
        nlon,nlat = (360,150)
        dlon,dlat = (1.0,1.0)
        extent = [-179.5,179.5,-59.5,89.5]

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
    hl,kl,ll = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE)

    # calculate Legendre polynomials
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # Setting units factor for cmwe: centimeters water equivalent
    # dfactor computes the degree dependent coefficients
    dfactor = gravtk.units(lmax=LMAX).harmonic(hl,kl,ll).cmwe
    # Gaussian smoothing
    if (RAD != 0):
        wt = 2.0*np.pi*gravtk.gauss_weights(RAD,LMAX)
    else:
        wt = np.ones((LMAX+1))

    # GRACE/GRACE-FO months of interest
    grace_month = sorted(set(np.arange(START_MON,END_MON+1)) - set(MISSING))
    # spatial object list for original and processed files
    original = []
    processed = []
    # for each GRACE/GRACE-FO month
    for t,gm in enumerate(grace_month):
        # calendar year and month
        calendar_year,calendar_month = gravtk.time.grace_to_calendar(gm)
        # GLDAS monthly data file from read_gldas_monthly.py
        args=(MODEL,SPACING,calendar_year,calendar_month,suffix[DATAFORM])
        f1='GLDAS_{0}{1}_TWC_{2:4d}_{3:02d}.{4}'.format(*args)
        # spherical harmonic data file
        args=(MODEL,SPACING,LMAX,order_str,gm,suffix[DATAFORM])
        f2='GLDAS_{0}{1}_TWC_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)
        # read data files for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            gldas_data = gravtk.spatial(spacing=[dlon,dlat],nlat=nlat,
                nlon=nlon,extent=extent).from_ascii(os.path.join(ddir,sub1,f1))
            gldas_Ylms = gravtk.harmonics().from_ascii(os.path.join(ddir,sub2,f2))
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            gldas_data = gravtk.spatial().from_netCDF4(os.path.join(ddir,sub1,f1))
            gldas_Ylms = gravtk.harmonics().from_netCDF4(os.path.join(ddir,sub2,f2))
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            gldas_data = gravtk.spatial().from_HDF5(os.path.join(ddir,sub1,f1))
            gldas_Ylms = gravtk.harmonics().from_HDF5(os.path.join(ddir,sub2,f2))
        # if destriping the monthly GLDAS spherical harmonic data
        if DESTRIPE:
            gldas_Ylms = gldas_Ylms.destripe()
        # convolve with degree dependent weighting
        gldas_Ylms.convolve(dfactor*wt)
        # convert spherical harmonics to output spatial grid
        gldas_grid = gldas_data.zeros_like()
        gldas_grid.time = np.copy(gldas_Ylms.time)
        gldas_grid.month = np.copy(gldas_Ylms.month)
        gldas_grid.data = gravtk.harmonic_summation(gldas_Ylms.clm, gldas_Ylms.slm,
            gldas_data.lon, gldas_data.lat, LMAX=LMAX, PLM=PLM).T
        # replace fill value points and certain vegetation types with 0
        gldas_data.replace_invalid(0.0, mask=combined_mask)
        gldas_grid.replace_invalid(0.0)
        gldas_grid.update_mask()
        # append original and processed data to list
        original.append(gldas_data)
        processed.append(gldas_grid)

    # create merged spatial objects from lists
    gldas_data = gravtk.spatial().from_list(original, clear=True)
    gldas_grid = gravtk.spatial().from_list(processed, clear=True)
    # calculate scaling factors and scaling factor errors
    gldas_kfactor = gldas_grid.kfactor(gldas_data)
    # calculate power of original data
    gldas_power = gldas_data.sum(power=2.0).power(0.5)

    # output file format
    file_format = 'GLDAS_{0}{1}_TWC_{2}_L{3:d}{4}{5}{6}_{7:03d}-{8:03d}.{9}'
    # attributes for output files
    attributes = {}
    attributes['title'] = gldas_products[MODEL]
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # Output scaling factor file
    f3 = file_format.format(MODEL,SPACING,'kfactor',LMAX,order_str,gw_str,
        ds_str,grace_month[0],grace_month[-1],suffix[DATAFORM])
    attributes['units'] = 'unitless'
    attributes['longname'] = 'Scaling_Factor'
    output_data(gldas_kfactor, FILENAME=os.path.join(ddir,sub2,f3),
        DATAFORM=DATAFORM, KEY='data', MODE=MODE, **attributes)
    # Output scaling factor error file
    f4 = file_format.format(MODEL,SPACING,'kf_error',LMAX,order_str,gw_str,
        ds_str,grace_month[0],grace_month[-1],suffix[DATAFORM])
    attributes['units'] = 'cmwe'
    attributes['longname'] = 'Scaling_Factor_Error'
    output_data(gldas_kfactor, FILENAME=os.path.join(ddir,sub2,f4),
        DATAFORM=DATAFORM, KEY='error', MODE=MODE, **attributes)
    # Output scaling factor power file
    f5 = file_format.format(MODEL,SPACING,'power',LMAX,order_str,gw_str,
        ds_str,grace_month[0],grace_month[-1],suffix[DATAFORM])
    attributes['units'] = 'cmwe'
    attributes['longname'] = 'Power'
    output_data(gldas_power, FILENAME=os.path.join(ddir,sub2,f5),
        DATAFORM=DATAFORM, KEY='data', MODE=MODE, **attributes)

# PURPOSE: wrapper function for outputting data to file
def output_data(data, FILENAME=None, KEY='data', DATAFORM=None,
    MODE=0o775, **kwargs):
    # copy data to output
    output = data.copy()
    setattr(output,'data',getattr(data,KEY))
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        output.to_ascii(FILENAME, date=False)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        output.to_netCDF4(FILENAME, date=False, **kwargs)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        output.to_HDF5(FILENAME, date=False, **kwargs)
    # change the permissions mode of the output file
    os.chmod(FILENAME, MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly GLDAS total water storage anomalies
            and monthly spherical harmonic coefficients to compute
            point scaling factors following Landerer and Swenson (2012)
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+', choices=gldas_products.keys(),
        help='GLDAS land surface model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # start and end years GRACE/GRACE-FO months
    parser.add_argument('--start','-s',
        type=int,default=4,
        help='Starting GRACE/GRACE-FO month')
    parser.add_argument('--end','-e',
        type=int,default=228,
        help='Ending GRACE/GRACE-FO month')
    # missing GRACE/GRACE-FO months
    MISSING = [6,7,18,109,114,125,130,135,140,141,146,151,
        156,162,166,167,172,177,178,182,187,188,189,190,191,
        192,193,194,195,196,197,200,201]
    parser.add_argument('--missing','-o',
        metavar='MISSING',type=int,nargs='+',default=MISSING,
        help='Missing GRACE/GRACE-FO months')
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
    # Use a decorrelation (destriping) filter
    parser.add_argument('--destripe','-d',
        default=False, action='store_true',
        help='Use decorrelation (destriping) filter')
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

    # for each GLDAS model
    for MODEL in args.model:
        # run program
        gldas_scaling_factors(args.directory, MODEL,
            args.start, args.end, args.missing,
            VERSION=args.version, SPACING=args.spacing, LMAX=args.lmax,
            MMAX=args.mmax, RAD=args.radius, DESTRIPE=args.destripe,
            LOVE_NUMBERS=args.love, REFERENCE=args.reference,
            DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
