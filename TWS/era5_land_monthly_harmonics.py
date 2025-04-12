#!/usr/bin/env python
u"""
era5_land_monthly_harmonics.py
Written by Tyler Sutterley (05/2023)

Reads monthly ERA5-Land total water storage anomalies and converts to
    spherical harmonic coefficients

snow_depth_water_equivalent (sd)
snow_cover (snowc): 0 -- 100 percent
skin_reservoir_content (src)
volumetric_soil_water_layer_1 (swvl1): 0 -- 7 cm
volumetric_soil_water_layer_2 (swvl2): 7 -- 28 cm
volumetric_soil_water_layer_3 (swvl3): 28 -- 100 cm
volumetric_soil_water_layer_4 (swvl4): 100 -- 289 cm

Time-averaged grid from a set yearly range subtracted from individual grids.

CALLING SEQUENCE:
    python era5_land_monthly_harmonics.py --lmax 60 --format netCDF4

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Years to run
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
    Written 04/2025
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

# PURPOSE: convert ERA5-Land terrestrial water storage data to spherical harmonics
def era5_land_monthly_harmonics(base_dir, YEARS,
        MASKS=None,
        LMAX=0,
        MMAX=None,
        LOVE_NUMBERS=0,
        REFERENCE=None,
        DATAFORM=None,
        MODE=0o775
    ):

    # directory models
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # ERA5-land products
    MODEL = 'ERA5-Land'
    ddir = base_dir.joinpath(MODEL)
    # Creating output subdirectory if it doesn't exist
    output_dir = ddir.joinpath(f'{MODEL}_TWS_CLM_L{LMAX:d}')
    output_dir.mkdir(mode=MODE, parents=True, exist_ok=True)
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # output dimensions
    nlat,nlon = (1801,3600)

    # attributes for output files
    attributes = {}
    attributes['product_version'] = MODEL
    attributes['product_name'] = 'TWS'
    attributes['product_type'] = 'gravity_field'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX
    # output string for both LMAX == MMAX and LMAX != MMAX cases
    order_str = 'M{MMAX:d}' if (MMAX != LMAX) else ''

    # ERA5-Land land/sea mask
    landmask_file = ddir.joinpath(f'lsm.nc')
    with netCDF4.Dataset(landmask_file, mode='r') as fileID:
        era5_land_mask = fileID.variables['lsm'][:].squeeze()
        glon = fileID.variables['longitude'][:].copy()
        glat = fileID.variables['latitude'][:].copy()
    # create mesh grid of latitude and longitude
    gridlon,gridlat = np.meshgrid(glon,glat)
    # create combined mask for invalid points
    combined_mask = (era5_land_mask <= 0.1)

    if MASKS:
        # read masks for reducing regions before converting to harmonics
        for mask_file in MASKS:
            logging.debug(str(mask_file))
            mask_file = pathlib.Path(mask_file).expanduser().absolute()
            with netCDF4.Dataset(mask_file, mode='r') as fileID:
                combined_mask |= fileID.variables['mask'][:].astype(bool)
    else:
        # use default masks for reducing regions before converting to harmonics
        # read Permafrost index file
        permafrost_file = ddir.joinpath(f'cpfrost.nc')
        logging.debug(str(permafrost_file))
        with netCDF4.Dataset(permafrost_file, mode='r') as fileID:
            fileID.set_auto_mask(False)
            permafrost_index = fileID.variables['pf'][0,:,:].copy()
        # 1: Continuous Permafrost
        # 2: Discontinuous Permafrost
        # 3: Isolated Permafrost
        # 4: Sporadic Permafrost
        # 5: Glaciated Area
        for invalid_keys in (1,5):
            combined_mask |= (permafrost_index == invalid_keys)
        # read Arctic mask file
        arctic_file = ddir.joinpath('carctic.nc')
        logging.debug(str(arctic_file))
        with netCDF4.Dataset(arctic_file, mode='r') as fileID:
            fileID.set_auto_mask(False)
            arctic_mask = fileID.variables['mask'][0,:,:].copy()
        # Arctic mask
        combined_mask |= (arctic_mask > 0)
        # read Glacier mask file
        glacier_file = ddir.joinpath(f'cicecap.nc')
        logging.debug(str(glacier_file))
        with netCDF4.Dataset(glacier_file, mode='r') as fileID:
            glacier_mask = fileID.variables['si10'][0,:,:].copy()
        # glacier mask
        combined_mask |= (glacier_mask > 0.0)

    # Earth Parameters
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

    # find input terrestrial water storage files
    regex_years = r'\d+' if (YEARS is None) else r'|'.join(map(str,YEARS))
    args = (MODEL, regex_years, suffix[DATAFORM])
    rx = re.compile(r'{0}-TWS-({1})-(\d+)\.{2}$'.format(*args))
    FILES = sorted([f for f in ddir.iterdir() if rx.match(f.name)])

    # for each input file
    for t,FILE in enumerate(FILES[:-1]):
        # extract year and month from file
        YY,MM = np.array(rx.findall(FILE.name).pop(), dtype=np.float64)

        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            M1 = gravtk.spatial().from_ascii(FILE, nlat=nlat, nlon=nlon)
            M2 = gravtk.spatial().from_ascii(FILES[t+1], nlat=nlat, nlon=nlon)
        elif (DATAFORM == 'netCDF4'):
            # netCDF4 (.nc)
            M1 = gravtk.spatial().from_netCDF4(FILE)
            M2 = gravtk.spatial().from_netCDF4(FILES[t+1])
        elif (DATAFORM == 'HDF5'):
            # HDF5 (.H5)
            M1 = gravtk.spatial().from_HDF5(FILE)
            M2 = gravtk.spatial().from_HDF5(FILES[t+1])
        # attributes for input files
        attributes['lineage'] = []
        attributes['lineage'].append(pathlib.Path(M1.filename).name)
        attributes['lineage'].append(pathlib.Path(M2.filename).name)

        # replace fill value points and certain vegetation types with 0
        M1.replace_invalid(0.0, mask=combined_mask)
        M2.replace_invalid(0.0, mask=combined_mask)
        # calculate 2-month moving average in cm w.e.
        # weighting by number of days in each month
        dpm = gravtk.time.calendar_days(int(YY))
        W = np.float64(dpm[(t+1) % 12] + dpm[t % 12])
        MASS = 100.0*(dpm[t % 12]*M1.data + dpm[(t+1) % 12]*M2.data)/W

        # convert to spherical harmonics
        era5_land_Ylms = gravtk.gen_stokes(MASS, glon, latitude_geocentric[:,0],
            LMAX=LMAX, MMAX=MMAX, PLM=PLM, LOVE=LOVE)
        # calculate date information
        era5_land_Ylms.time, = gravtk.time.convert_calendar_decimal(YY,MM)
        # calculate GRACE/GRACE-FO month
        era5_land_Ylms.month = gravtk.time.calendar_to_grace(YY,MM)
        # add attributes to output harmonics
        era5_land_Ylms.attributes['ROOT'] = attributes

        # output spherical harmonic data file
        args=(MODEL, LMAX, order_str, era5_land_Ylms.month, suffix[DATAFORM])
        FILE = '{0}_TWS_CLM_L{1:d}{2}_{3:03d}.{4}'.format(*args)
        output_file = output_dir.joinpath(FILE)
        era5_land_Ylms.to_file(output_file, format=DATAFORM)
        # change the permissions mode of the output file to MODE
        output_file.chmod(mode=MODE)

    # Output date ascii file
    output_date_file = output_dir.joinpath(f'{MODEL}_TWS_DATES.txt')
    fid1 = output_date_file.open(mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = output_dir.joinpath('index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # find all available output files
    args = (MODEL, LMAX, order_str, suffix[DATAFORM])
    rx = re.compile(r'{0}_TWS_CLM_L{1:d}{2}_([-]?\d+).{3}'.format(*args))
    # find all output harmonic files (not just ones created in run)
    output_files = [fi for fi in output_dir.iterdir() if rx.match(fi.name)]
    for fi in sorted(output_files):
        # extract GRACE month
        grace_month, = np.array(rx.findall(fi.name), dtype=int)
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
        description="""Reads monthly ERA5-Land total water storage anomalies
            and converts to spherical harmonic coefficients
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
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

    # run program
    era5_land_monthly_harmonics(args.directory, args.year,
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
