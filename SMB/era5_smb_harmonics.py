#!/usr/bin/env python
u"""
era5_smb_harmonics.py
Written by Tyler Sutterley (03/2023)
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
import re
import logging
import netCDF4
import pathlib
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
    ddir = pathlib.Path(ddir).expanduser().absolute()
    d1 = ddir.joinpath('ERA5-Cumul-P-E-{0:4d}-{1:4d}'.format(*RANGE))
    # Creating output subdirectory if it doesn't exist
    prefix = f'{REGION}_' if REGION else ''
    d2 = ddir.joinpath(f'{prefix}ERA5_CUMUL_P-E_CLM_L{LMAX:d}')
    d2.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output data file format and title
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # attributes for output files
    attributes = {}
    attributes['project'] = 'ECMWF atmospheric reanalysis'
    attributes['title'] = 'ERA5 Precipitation minus Evaporation'
    attributes['product_name'] = 'P-E'
    attributes['source'] = ', '.join(['tp','e'])
    attributes['product_type'] = 'gravity_field'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

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
        logging.debug(str(mask_file))
        mask_file = pathlib.Path(mask_file).expanduser().absolute()
        fileID = netCDF4.Dataset(mask_file, mode='r')
        input_mask |= fileID.variables['mask'][:].astype(bool)
        fileID.close()

    # get reference parameters for ellipsoid
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [m]
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

    # find input files from era5_smb_cumulative.py
    regex_years = r'\d{4}' if (YEARS is None) else '|'.join(map(str,YEARS))
    rx = re.compile(r'ERA5\-Cumul\-P-E\-({0})\.nc$'.format(regex_years))
    input_files = sorted([f for f in d1.iterdir() if rx.match(f.name)])

    # create list of yearly ERA5 files
    spatial_list = []
    # for each input file
    for t,input_file in enumerate(input_files):
        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            era5_data = gravtk.spatial().from_ascii(input_file,
                spacing=[dlon,dlat], nlat=nlat, nlon=nlon,
                extent=extent)
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            era5_data = gravtk.spatial().from_netCDF4(input_file, varname='SMB')
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            era5_data = gravtk.spatial().from_HDF5(input_file, varname='SMB')
        # if reducing to a region of interest before converting to harmonics
        if np.any(input_mask):
            # replace fill value points and masked points with 0
            era5_data.replace_invalid(0.0, mask=np.logical_not(input_mask))
        else:
            # replace fill value points points with 0
            era5_data.replace_invalid(0.0)
        # extend monthly values to spatial list
        spatial_list.extend(era5_data)
    # convert to combined data cube and clear memory from spatial list
    era5_data = gravtk.spatial().from_list(spatial_list, clear=True)
    nlat,nlon,nt = era5_data.shape

    # for each month of data
    for i in range(nt-1):
        # convert data to mm w.e.
        M1 = era5_data.index(i).scale(1000.0)
        M2 = era5_data.index(i+1).scale(1000.0)
        # attributes for input files
        attributes['lineage'] = []
        attributes['lineage'].append(pathlib.Path(M1.filename).name)
        attributes['lineage'].append(pathlib.Path(M2.filename).name)
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
        # add attributes to output harmonics
        era5_Ylms.attributes['ROOT'] = attributes
        # output spherical harmonic data file
        args = (LMAX, order_str, era5_Ylms.month, suffix[DATAFORM])
        FILE = 'ERA5_CUMUL_P-E_CLM_L{0:d}{1}_{2:03d}.{3}'.format(*args)
        output_file = d2.joinpath(FILE)
        era5_Ylms.to_file(output_file, format=DATAFORM)
        # change the permissions mode of the output file to MODE
        output_file.chmod(mode=MODE)

    # Output date ascii file
    output_date_file = d2.joinpath('ERA5_SMB_DATES.txt')
    fid1 = output_date_file.open(mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = d2.joinpath('index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # find all available output files
    args = (LMAX, order_str, suffix[DATAFORM])
    output_pattern = r'ERA5_CUMUL_P-E_CLM_L{0:d}{1}_([-]?\d+).{2}'
    output_regex = re.compile(output_pattern.format(*args), re.VERBOSE)
    # find all output harmonic files (not just ones created in run)
    output_files = [f for f in d2.iterdir() if re.match(output_regex,f.name)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi.name), dtype=int)
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
        description="""Reads monthly ERA5 surface mass balance
            anomalies and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
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

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program with parameters
    era5_smb_harmonics(args.directory, args.year, RANGE=args.mean,
        REGION=args.region, MASKS=args.mask, LMAX=args.lmax, MMAX=args.mmax,
        LOVE_NUMBERS=args.love, REFERENCE=args.reference,
        DATAFORM=args.format, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
