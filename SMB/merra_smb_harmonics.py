#!/usr/bin/env python
u"""
merra_smb_harmonics.py
Written by Tyler Sutterley (03/2023)
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
    -n X, --love X: Treatment of the Love Love numbers
        0: Han and Wahr (1995) values from PREM
        1: Gegout (2005) values from PREM
        2: Wang et al. (2012) values from PREM
        3: Wang et al. (2012) values from PREM with hard sediment
        4: Wang et al. (2012) values from PREM with soft sediment
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
    datum.py: calculate reference parameters for common ellipsoids
    gen_stokes.py: converts a spatial field into a series of spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
        updated inputs to spatial from_ascii function
        use spatial function for calculating geocentric latitude
    Updated 02/2023: use love numbers class with additional attributes
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
import re
import copy
import logging
import netCDF4
import pathlib
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
    VERSION = '5.12.4'
    cumul_sub = f'{PRODUCT}.{VERSION}.CUMUL.{RANGE[0]:d}.{RANGE[1]:d}'
    input_dir = ddir.joinpath(cumul_sub)
    # Creating output subdirectory if it doesn't exist
    prefix = f'{REGION}_' if REGION else ''
    output_sub = f'{prefix}{PRODUCT}_{VERSION}_CUMUL_CLM_L{LMAX:d}'
    output_dir = ddir.joinpath(output_sub)
    output_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
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
    attributes['institution'] = 'NASA Goddard Space Flight Center (GSFC)'
    attributes['project'] = 'MERRA2'
    attributes['title'] = copy.copy(merra_products[PRODUCT])
    attributes['source'] = ', '.join(merra_sources[PRODUCT])
    attributes['product_name'] = PRODUCT
    attributes['product_version'] = VERSION
    attributes['product_type'] = 'gravity_field'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

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
        mask_file = pathlib.Path(mask_file).expanduser().absolute()
        fileID = netCDF4.Dataset(mask_file, mode='r')
        input_mask |= fileID.variables['mask'][:].astype(bool)
        fileID.close()

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [cm]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = mdlhmc.spatial.geocentric_latitude(gridlon, gridlat,
        a_axis=a_axis, flat=flat)
    # colatitude in radians
    theta = (90.0 - latitude_geocentric[:,0])*np.pi/180.0

    # read load love numbers
    LOVE = gravtk.load_love_numbers(LMAX, LOVE_NUMBERS=LOVE_NUMBERS,
        REFERENCE=REFERENCE, FORMAT='class')
    # add attributes for earth parameters
    attributes['earth_model'] = LOVE.model
    attributes['earth_love_numbers'] = LOVE.citation
    attributes['reference_frame'] = LOVE.reference
    # add attributes for maximum degree and order
    attributes['max_degree'] = LMAX
    attributes['max_order'] = MMAX

    # calculate Legendre polynomials
    PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(theta))

    # find input files from merra_smb_cumulative.py
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (PRODUCT, regex_years, suffix[DATAFORM])
    regex_pattern = r'MERRA2_(\d+).tavgM_2d_{0}_cumul_Nx.(({1})(\d{{2}})).{2}$'
    rx = re.compile(regex_pattern.format(*args), re.VERBOSE)
    # will be out of order for 2020 due to September reprocessing
    FILES = sorted([f.name for f in input_dir.iterdir() if rx.match(f.name)])
    # sort files by month
    indices = np.argsort([rx.match(f1).group(2) for f1 in FILES])
    FILES = [FILES[indice] for indice in indices]
    # remove files that needed to be reprocessed
    INVALID = []
    INVALID.append('MERRA2_400.tavgM_2d_{0}_cumul_Nx.202009.{2}'.format(*args))
    INVALID.append('MERRA2_400.tavgM_2d_{0}_cumul_Nx.202106.{2}'.format(*args))
    INVALID.append('MERRA2_400.tavgM_2d_{0}_cumul_Nx.202107.{2}'.format(*args))
    INVALID.append('MERRA2_400.tavgM_2d_{0}_cumul_Nx.202109.{2}'.format(*args))
    if (set(INVALID) & set(FILES)):
        logging.warning("Reprocessed file found in list")
        FILES = sorted(set(FILES) - set(INVALID))

    # for each input file
    for t,fi in enumerate(FILES[:-1]):
        # extract parameters from input flux file
        MOD,_,YY,MM = rx.findall(fi).pop()
        f1 = input_dir.joinpath(fi)
        f2 = input_dir.joinpath(FILES[t+1])
        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            M1 = gravtk.spatial().from_ascii(f1, spacing=[dlon,dlat],
                nlat=nlat, nlon=nlon, extent=extent)
            M2 = gravtk.spatial().from_ascii(f2, spacing=[dlon,dlat],
                nlat=nlat, nlon=nlon, extent=extent)
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            M1 = gravtk.spatial().from_netCDF4(f1, varname=PRODUCT)
            M2 = gravtk.spatial().from_netCDF4(f2, varname=PRODUCT)
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            M1 = gravtk.spatial().from_HDF5(f1, varname=PRODUCT)
            M2 = gravtk.spatial().from_HDF5(f2, varname=PRODUCT)
        # attributes for input files
        attributes['lineage'] = []
        attributes['lineage'].append(f1.name)
        attributes['lineage'].append(f2.name)

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
        # add attributes to output harmonics
        merra_Ylms.attributes['ROOT'] = attributes
        # output spherical harmonic data file
        args = (MOD,PRODUCT,LMAX,order_str,merra_Ylms.month,suffix[DATAFORM])
        FILE='MERRA2_{0}_tavgM_2d_{1}_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)
        output_file = output_dir.joinpath(FILE)
        merra_Ylms.to_file(output_file, format=DATAFORM)
        # change the permissions mode of the output file to MODE
        output_file.chmod(mode=MODE)

    # Output date ascii file
    output_date_file = output_dir.joinpath(f'MERRA2_{PRODUCT}_DATES.txt')
    fid1 = output_date_file.open(mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = output_dir.joinpath('index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # find all available output files
    args = (PRODUCT, LMAX, order_str, suffix[DATAFORM])
    output_pattern = r'MERRA2_(\d+)_tavgM_2d_{0}_CLM_L{1:d}{2}_([-]?\d+).{3}'
    output_regex = re.compile(output_pattern.format(*args), re.VERBOSE)
    # find all output harmonic files (not just ones created in run)
    output_files = [f for f in output_dir.iterdir()
        if re.match(output_regex,f.name)]
    for fi in sorted(output_files):
        # extract GRACE month
        MOD,grace_month = np.array(re.findall(output_regex,fi.name).pop(), dtype=int)
        YY,MM = gravtk.time.grace_to_calendar(grace_month)
        tdec, = gravtk.time.convert_calendar_decimal(YY, MM)
        # print date, GRACE month and calendar month to date file
        fid1.write(f'{tdec:11.6f} {grace_month:03d} {MM:02.0f}\n')
        # print output file to index
        full_output_file = gravtk.spatial().compressuser(fi)
        print(full_output_file, file=fid2)
    # close the date and index files
    fid1.close()
    fid2.close()
    # set the permissions level of the output date and index files to MODE
    output_date_file.chmod(mode=MODE)
    output_index_file.chmod(mode=MODE)

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
        type=pathlib.Path, default=pathlib.Path.cwd(),
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
        type=pathlib.Path,
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
    # 3: Wang et al. (2012) values from PREM with hard sediment
    # 4: Wang et al. (2012) values from PREM with soft sediment
    parser.add_argument('--love','-n',
        type=int, default=0, choices=[0,1,2,3,4],
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
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
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
