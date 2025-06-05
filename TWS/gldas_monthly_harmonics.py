#!/usr/bin/env python
u"""
gldas_monthly_harmonics.py
Written by Tyler Sutterley (04/2025)

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
    Updated 04/2025: don't use auto mask with netCDF4 files (set_auto_mask)
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: add root attributes to output netCDF4 and HDF5 files
        updated inputs to spatial from_ascii function
        use spatial function for calculating geocentric latitude
    Updated 02/2023: use love numbers class with additional attributes
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
import re
import logging
import netCDF4
import pathlib
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
def gldas_monthly_harmonics(base_dir, MODEL, YEARS,
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
    # directory for GLDAS models
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # subdirectory for model monthly products at spacing for version
    d1 = base_dir.joinpath(f'GLDAS_{MODEL}{SPACING}_{TEMPORAL}{V2}')
    # Creating output subdirectory if it doesn't exist
    d2 = base_dir.joinpath(f'GLDAS_{MODEL}{SPACING}{V1}_TWC_CLM_L{LMAX:d}')
    d2.mkdir(mode=MODE, parents=True, exist_ok=True)

    # attributes for output files
    attributes = {}
    attributes['institution'] = 'NASA Goddard Space Flight Center (GSFC)'
    attributes['project'] = 'Global Land Data Assimilation System (GLDAS)'
    attributes['product_version'] = f'{MODEL} v{VERSION}'
    attributes['product_name'] = 'TWC'
    attributes['product_type'] = 'gravity_field'
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

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
    landmask_file = base_dir.joinpath(f'GLDASp5_landmask_{SPACING}d.nc4')
    with netCDF4.Dataset(landmask_file, mode='r') as fileID:
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
            logging.debug(str(mask_file))
            mask_file = pathlib.Path(mask_file).expanduser().absolute()
            with netCDF4.Dataset(mask_file, mode='r') as fileID:
                combined_mask |= fileID.variables['mask'][:].astype(bool)
    else:
        # use default masks for reducing regions before converting to harmonics
        # mask combining vegetation index, permafrost index and Arctic mask
        # read vegetation index file
        vegetation_file = base_dir.joinpath(f'modmodis_domveg20_{SPACING}.nc')
        logging.debug(str(vegetation_file))
        with netCDF4.Dataset(vegetation_file, mode='r') as fileID:
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
        permafrost_file = base_dir.joinpath(f'permafrost_mod44w_{SPACING}.nc')
        logging.debug(str(permafrost_file))
        with netCDF4.Dataset(permafrost_file, mode='r') as fileID:
            fileID.set_auto_mask(False)
            permafrost_index = fileID.variables['mask'][:]
        # 1: Continuous Permafrost
        # 2: Discontinuous Permafrost
        # 3: Isolated Permafrost
        # 4: Sporadic Permafrost
        # 5: Glaciated Area
        for invalid_keys in (1,5):
            combined_mask |= (permafrost_index == invalid_keys)
        # read Arctic mask file
        arctic_file = base_dir.joinpath(f'arcticmask_mod44w_{SPACING}.nc')
        logging.debug(str(arctic_file))
        with netCDF4.Dataset(arctic_file, mode='r') as fileID:
            fileID.set_auto_mask(False)
            arctic_mask = fileID.variables['mask'][:].astype(bool)
        # arctic mask
        combined_mask |= arctic_mask[:,:]

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
    args = (MODEL, SPACING, regex_years, suffix[DATAFORM])
    rx = re.compile(r'GLDAS_{0}{1}_TWC_({2})_(\d+)\.{3}$'.format(*args))
    FILES = sorted([f for f in d1.iterdir() if rx.match(f.name)])

    # for each input file
    for t,FILE in enumerate(FILES[:-1]):
        # extract year and month from file
        YY,MM = np.array(rx.findall(FILE.name).pop(), dtype=np.float64)

        # read data file for data format
        if (DATAFORM == 'ascii'):
            # ascii (.txt)
            M1 = gravtk.spatial().from_ascii(FILE,
                spacing=[dlon,dlat], nlat=nlat, nlon=nlon, extent=extent)
            M2 = gravtk.spatial().from_ascii(FILES[t+1],
                spacing=[dlon,dlat], nlat=nlat, nlon=nlon, extent=extent)
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
        # add attributes to output harmonics
        gldas_Ylms.attributes['ROOT'] = attributes

        # output spherical harmonic data file
        args=(MODEL,SPACING,LMAX,order_str,gldas_Ylms.month,suffix[DATAFORM])
        FILE = 'GLDAS_{0}{1}_TWC_CLM_L{2:d}{3}_{4:03d}.{5}'.format(*args)
        output_file = d2.joinpath(FILE)
        gldas_Ylms.to_file(output_file, format=DATAFORM)
        # change the permissions mode of the output file to MODE
        output_file.chmod(mode=MODE)

    # Output date ascii file
    output_date_file = d2.joinpath(f'GLDAS_{MODEL}{SPACING}_TWC_DATES.txt')
    fid1 = output_date_file.open(mode='w', encoding='utf8')
    # date file header information
    print('{0:8} {1:^6} {2:^5}'.format('Mid-date','GRACE','Month'), file=fid1)
    # index file listing all output spherical harmonic files
    output_index_file = d2.joinpath('index.txt')
    fid2 = output_index_file.open(mode='w', encoding='utf8')
    # find all available output files
    args = (MODEL, SPACING, LMAX, order_str, suffix[DATAFORM])
    output_regex=r'GLDAS_{0}{1}_TWC_CLM_L{2:d}{3}_([-]?\d+).{4}'.format(*args)
    # find all output harmonic files (not just ones created in run)
    output_files = [fi for fi in d2.iterdir() if re.match(output_regex,fi.name)]
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
        type=pathlib.Path, default=pathlib.Path.cwd(),
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
