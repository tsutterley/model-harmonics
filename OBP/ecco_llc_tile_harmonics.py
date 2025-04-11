#!/usr/bin/env python
u"""
ecco_llc_tile_harmonics.py
Written by Tyler Sutterley (05/2023)
Reads monthly ECCO ocean bottom pressure anomalies from LLC tiles
    and converts to spherical harmonic coefficients

ECCO Version 4, Revision 4
    https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_monthly

ECCO Version 5, Alpha
    https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly

INPUTS:
    ECCO version 4 or 5 models
        V4r4: Version 4, Revision 4
        V5alpha: Version 5, Alpha release

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to run
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
    -F X, --format X: Output data format
        ascii
        netcdf
        HDF5
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    read_love_numbers.py: reads Load Love Numbers from Han and Wahr (1995)
    datum.py: calculate reference parameters for common ellipsoids
    ref_ellipsoid.py: calculate reference parameters for common ellipsoids
    norm_gravity.py: calculates the normal gravity for locations on an ellipsoid
    gen_point_pressure.py: converts pressure point values to spherical harmonics
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: add attributes to output netCDF4 and HDF5 files
    Updated 02/2023: use love numbers class with additional attributes
    Updated 12/2022: single implicit import of spherical harmonic tools
        use constants class in place of geoid-toolkit ref_ellipsoid
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: use wrapper function for reading load Love numbers
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
        use output harmonic file wrapper routine to write to file
    Updated 09/2021: use GRACE/GRACE-FO month to calendar month converters
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Written 02/2021
"""
from __future__ import print_function

import sys
import re
import logging
import pathlib
import netCDF4
import argparse
import datetime
import numpy as np
import gravity_toolkit as gravtk
import geoid_toolkit as geoidtk
import model_harmonics as mdlhmc

# PURPOSE: convert monthly ECCO OBP data to spherical harmonics
def ecco_llc_tile_harmonics(ddir, MODEL, YEARS, LMAX=0, MMAX=None,
    LOVE_NUMBERS=0, REFERENCE=None, DATAFORM=None, MODE=0o775):

    # input and output subdirectories
    ddir = pathlib.Path(ddir).expanduser().absolute()
    d1 = ddir.joinpath(f'ECCO_{MODEL}_AveRmvd_OBP','nctiles_monthly')
    d2 = ddir.joinpath(f'ECCO_{MODEL}_AveRmvd_OBP_CLM_L{LMAX:d}')
    # Creating subdirectory if it doesn't exist
    d2.mkdir(mode=MODE, parents=True, exist_ok=True)

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''
    # output file format
    output_file_format = 'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}{2}_{3:03d}.{4}'
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # input variable names for each model
    if (MODEL == 'V4r4'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        AREANAME = 'rA'
        Nt,Nj,Ni = (13,90,90)
    elif (MODEL == 'V5alpha'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        AREANAME = 'rA'
        Nt,Nj,Ni = (13,270,270)

    # read ECCO tile grid file
    input_invariant_file = ddir.joinpath(f'ECCO-{MODEL}',
        'nctiles_monthly', 'ECCO-GRID.nc')
    invariant = ncdf_invariant(input_invariant_file,
        lon=LONNAME, lat=LATNAME, depth=ZNAME, area=AREANAME)
    # read ECCO tile geoid height file from ecco_geoid_llc_tiles.py
    input_geoid_file = ddir.joinpath(f'ECCO-{MODEL}',
        'nctiles_monthly', 'ECCO-EGM2008.nc')
    geoid_undulation = ncdf_geoid(input_geoid_file)
    # use geoid and depth to calculate bathymetry
    bathymetry = geoid_undulation - invariant['depth']
    # attributes for output files
    attributes = {}
    attributes['institution'] = 'NASA Jet Propulsion Laboratory (JPL)'
    PROJECT = 'Estimating the Circulation and Climate of the Ocean (ECCO)'
    attributes['project'] = PROJECT
    attributes['product_version'] = MODEL
    attributes['product_name'] = VARNAME
    attributes['product_type'] = 'gravity_field'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # Earth Parameters
    ellipsoid_params = mdlhmc.datum(ellipsoid='WGS84')
    # semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude and longitude in radians
    latitude_geodetic_rad = np.pi*invariant['lat']/180.0
    longitude_rad = np.pi*invariant['lon']/180.0
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.*np.sin(latitude_geodetic_rad)**2.)
    # calculate X, Y and Z from geodetic latitude and longitude
    X = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.cos(longitude_rad)
    Y = (N+bathymetry)*np.cos(latitude_geodetic_rad)*np.sin(longitude_rad)
    Z = (N * (1.0 - ecc1**2.0) + bathymetry) * np.sin(latitude_geodetic_rad)
    R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = 180.0*np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/np.pi
    # colatitude in radians
    theta = (90.0 - latitude_geocentric)*np.pi/180.0

    # calculate normal gravity at latitudes and bathymetry
    gamma_h,dgamma_dh = geoidtk.norm_gravity(latitude_geocentric,
        bathymetry, 'WGS84')

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

    # regular expression pattern to find files and extract dates
    regex_years = r'\d+' if (YEARS is None) else '|'.join(map(str,YEARS))
    args = (MODEL, regex_years, suffix[DATAFORM])
    rx = re.compile(r'ECCO_{0}_AveRmvd_OBP_({1})_(\d+).{2}$'.format(*args))

    # find input ECCO OBP files
    FILES = [fi for fi in d1.iterdir() if rx.match(fi.name)]
    # for each input file
    for f in sorted(FILES):
        # extract dates from file
        year,month = np.array(rx.findall(f.name).pop(), dtype=np.int64)
        # read input data file
        obp_data = {}
        with netCDF4.Dataset(f, mode='r') as fileID:
            for key,val in fileID.variables.items():
                obp_data[key] = val[:].copy()
        # allocate for output spherical harmonics
        obp_Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
        obp_Ylms.clm = np.zeros((LMAX+1,MMAX+1))
        obp_Ylms.slm = np.zeros((LMAX+1,MMAX+1))
        # copy date information
        obp_Ylms.time = np.copy(obp_data['time'])
        obp_Ylms.month = gravtk.time.calendar_to_grace(year,month)
        # calculate harmonics for each spatial tile
        # add to the total set of harmonics
        for k in range(Nt):
            # valid indices for for tile
            indj,indi = np.nonzero(~obp_data[VARNAME].mask[k,:,:])
            # reduce to valid
            obp_tile = obp_data[VARNAME].data[k,indj,indi]
            gamma_tile = gamma_h[k,indj,indi]
            rad_tile = R[k,indi,indj]
            # grid point areas (m^2)
            area = invariant['area'][k,indj,indi]
            # geocentric latitude and longitude
            lat = latitude_geocentric[k,indj,indi]
            lon = invariant['lon'][k,indj,indi]
            # calculate harmonics from pressure/gravity ratio
            Ylms = mdlhmc.gen_point_pressure(obp_tile, gamma_tile, rad_tile,
                lon, lat, AREA=area, LMAX=LMAX, MMAX=MMAX, LOVE=LOVE)
            # add tile Ylms to total for month
            obp_Ylms.add(Ylms)
        # attributes for input files
        attributes['lineage'] = []
        attributes['lineage'].append(input_invariant_file.name)
        attributes['lineage'].append(input_geoid_file.name)
        attributes['lineage'].append(f.name)
        # add attributes to output harmonics
        obp_Ylms.attributes['ROOT'] = attributes
        # output spherical harmonic data file
        args = (MODEL, LMAX, order_str, obp_Ylms.month, suffix[DATAFORM])
        output_file = d2.joinpath(output_file_format.format(*args))
        obp_Ylms.to_file(output_file,
            format=DATAFORM, date=True)
        # change the permissions mode of the output file to MODE
        output_file.chmod(mode=MODE)

    # Output date ascii file
    output_date_file = d2.joinpath(f'ECCO_{MODEL}_OBP_DATES.txt')
    fid1 = output_date_file.open( mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = d2.joinpath('index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # find all available output files
    args = (MODEL, LMAX, suffix[DATAFORM])
    output_regex=r'ECCO_{0}_AveRmvd_OBP_CLM_L{1:d}_([-]?\d+).{2}'.format(*args)
    # find all output harmonic files (not just ones created in run)
    output_files = [fi for fi in d2.iterdir() if re.match(output_regex,fi.name)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(re.findall(output_regex,fi.name),dtype=np.int64)
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

# PURPOSE: read ECCO invariant grid file
def ncdf_invariant(invariant_file, **kwargs):
    # output dictionary with invariant parameters
    invariant = {}
    # open netCDF4 file for reading
    with netCDF4.Dataset(invariant_file, mode='r') as fileID:
        # extract latitude, longitude, depth, area and valid mask
        for key,val in kwargs.items():
            invariant[key] = fileID.variables[val][:].copy()
    # return the invariant parameters
    return invariant

# PURPOSE: read geoid height netCDF4 files
def ncdf_geoid(FILENAME):
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        geoid_undulation = fileID.variables['geoid'][:,:,:].copy()
    return geoid_undulation

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO ocean bottom pressure
            anomalies from LLC tiles and converts to spherical harmonic
            coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['V5alpha'],
        choices=['V4r4','V5alpha'],
        help='ECCO Version 4 or 5 Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
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
    # Output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Output data format')
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

    # for each ECCO model
    for MODEL in args.model:
        # run program
        ecco_llc_tile_harmonics(args.directory, MODEL, args.year,
            LMAX=args.lmax, MMAX=args.mmax, LOVE_NUMBERS=args.love,
            REFERENCE=args.reference, DATAFORM=args.format,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
