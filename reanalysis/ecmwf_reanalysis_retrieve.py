#!/usr/bin/env python
"""
ecmwf_reanalysis_retrieve.py (07/2026)
Retrieves ERA5 reanalysis netCDF4 datasets from the ECMWF datastores API

levtype: surface (sfc), pressure levels (levtype:pl), model levels (ml)
type: an (analysis), fc (forecast)
Available parameters: http://apps.ecmwf.int/codes/grib/param-db

Get UID and API key from the ECMWF portal
https://ecds.ecmwf.int/how-to-api
https://ecds.ecmwf.int/profile
Write credentials into the .ecmwfdatastoresrc configuration file
Alternatively write as ECMWF_DATASTORES_URL and ECMWF_DATASTORES_KEY
    environmental variables

COMMAND LINE OPTIONS:
    -U X, --api-url: ECMWF datastores api url
    -K X, --api-key: ECMWF datastores api key
    -D X, --directory X: Working data directory
    -Y X, --year X: Year to retrieve
    -S X, --surface X: Retrieve model surface variables
        MSL: mean sea level pressure field
        SP: surface pressure field
        T2m: 2-metre temperature field
        P-E: precipitation and evaporation fields
    -L, --level: Retrieve the model level variables
    -I, --invariant: Retrieve the model invariant parameters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    ecmwf-datastores-client: ECMWF Data Stores Service API Python client
        https://ecmwf.github.io/ecmwf-datastores-client/

UPDATE HISTORY:
    Updated 07/2026: use new ecmwf.datastores API for retrieving data
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 07/2021: added option for retrieving the model level variables
    Updated 03/2021: added mean sea level pressure (msl) field as output
        automatically update years to run based on current time
    Updated 01/2021: added command line options for ECMWF api credentials
    Updated 12/2020: using argparse to set parameters
    Updated 07/2018: close the server connection after completion of program
        added 2-metre temperature (T2m) field as output
    Updated 03/2018: added option --mode to set the permissions level
    Written 03/2018
"""

from __future__ import print_function

import os
import time
import pathlib
import argparse
from ecmwf.datastores import Client


# PURPOSE: retrieve ECMWF data using the ECMWF datastores API
def ecmwf_reanalysis_retrieve(
    base_dir,
    client,
    YEAR,
    SURFACE=[],
    LEVEL=False,
    INVARIANT=True,
    MODE=0o775,
):
    # verify input data directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()

    # parameters for ERA5 dataset
    MODEL = 'ERA5'
    model_class = 'ea'
    model_dataset = 'era5'
    model_grid = '0.25/0.25'
    # surface variables
    surface_variable_dict = {}
    # mean sea level pressure field
    surface_variable_dict['MSL'] = 'mean_sea_level_pressure'
    # surface pressure field
    surface_variable_dict['SP'] = 'surface_pressure'
    # 2-metre temperature field
    surface_variable_dict['T2m'] = '2m_temperature'
    # precipitation and evaporation fields
    surface_variable_dict['P-E'] = [
        'total_precipitation',
        'convective_precipitation',
        'large_scale_precipitation',
        'snowfall',
        'evaporation',
    ]
    # model levels
    model_levelist = '/'.join([f'{l:d}' for l in range(1, 137 + 1)])

    # setup output directory and recursively create if currently non-existent
    ddir = base_dir.joinpath(MODEL)
    ddir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # for each year
    for y in YEAR:
        # months to retrieve
        months = [f'{m + 1:02d}' for m in range(12)]
        # monthly dates to retrieve
        d = '/'.join([f'{y:4d}{m}{1:02d}' for m in months])

        # for each surface variable to retrieve
        for surf in SURFACE:
            output_file = ddir.joinpath(f'{MODEL}-Monthly-{surf}-{y:4d}.nc')
            collection_id = 'reanalysis-era5-single-levels-monthly-means'
            request = {
                'year': str(y),
                'month': months,
                'time': '00:00',
                'grid': model_grid,
                'variable': surface_variable_dict[surf],
                'format': 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
            }
            # submit request to ECMWF datastores API
            remote = client.submit(collection_id, request)
            remote.download(str(output_file))
            # change the permissions mode to MODE
            output_file.chmod(mode=MODE)

        # if retrieving the model level data
        if LEVEL:
            # retrieve model temperature and specific humidity
            output_file = ddir.joinpath(f'{MODEL}-Monthly-Levels-{y:4d}.nc')
            collection_id = 'reanalysis-era5-complete'
            request = {
                'class': model_class,
                'dataset': model_dataset,
                'date': d,
                'expver': '1',
                'grid': model_grid,
                'levelist': model_levelist,
                'levtype': 'ml',
                'param': '130.128/133.128',
                'stream': 'moda',
                'type': 'an',
                'format': 'netcdf',
            }
            # submit request to ECMWF datastores API
            remote = client.submit(collection_id, request)
            remote.download(str(output_file))
            # change the permissions mode to MODE
            output_file.chmod(mode=MODE)

    # if retrieving the model invariant parameters
    if INVARIANT:
        output_file = ddir.joinpath(f'{MODEL}-Invariant-Parameters.nc')
        collection_id = 'reanalysis-era5-single-levels-monthly-means'
        request = {
            'class': model_class,
            'dataset': model_dataset,
            'year': '1979',
            'month': '01',
            'time': '00:00',
            'expver': '1',
            'grid': model_grid,
            'levtype': 'sfc',
            'variable': [
                'angle_of_sub_gridscale_orography',
                'anisotropy_of_sub_gridscale_orography',
                'geopotential',
                'land_sea_mask',
                'orography',
                'slope_of_sub_gridscale_orography',
                'standard_deviation_of_filtered_subgrid_orography',
                'standard_deviation_of_orography',
            ],
            'product_type': 'monthly_averaged_reanalysis',
            'format': 'netcdf',
        }
        # submit request to ECMWF datastores API
        remote = client.submit(collection_id, request)
        remote.download(str(output_file))
        # change the permissions mode to MODE
        output_file.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Retrieves ERA5 reanalysis netCDF4 datasets
            from the ECMWF datastores API
            """
    )
    # command line parameters
    # ECMWF api credentials
    parser.add_argument(
        '--api-url',
        '-U',
        type=str,
        default=os.environ.get('ECMWF_DATASTORES_URL'),
        help='ECMWF datastores api url',
    )
    parser.add_argument(
        '--api-key',
        '-K',
        type=str,
        default=os.environ.get('ECMWF_DATASTORES_KEY'),
        help='ECMWF datastores api key',
    )
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # years to retrieve
    now = time.gmtime()
    parser.add_argument(
        '--year',
        '-Y',
        type=int,
        nargs='+',
        default=range(2000, now.tm_year + 1),
        help='Model years to retrieve',
    )
    # retrieve model surface variables
    # MSL: mean sea level pressure field
    # SP: surface pressure field
    # T2m: 2-metre temperature field
    # P-E: Precipitation and Evaporation fields
    choices = ['MSL', 'SP', 'T2m', 'P-E']
    parser.add_argument(
        '--surface',
        '-S',
        type=str,
        nargs='?',
        metavar='VARIABLE',
        choices=choices,
        const=[],
        default=['SP'],
        help='Retrieve model surface variables',
    )
    # retrieve the model level variables
    parser.add_argument(
        '--level',
        '-L',
        default=False,
        action='store_true',
        help='Retrieve model level variables',
    )
    # retrieve the model invariant parameters
    parser.add_argument(
        '--invariant',
        '-I',
        default=False,
        action='store_true',
        help='Retrieve model invariant parameters',
    )
    # connection timeout
    parser.add_argument(
        '--timeout',
        '-t',
        type=int,
        default=360,
        help='Timeout in seconds for blocking operations',
    )
    # permissions mode of the directories and files retrieved
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files retrieved',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # open connection with ECMWF datastores api server
    client = Client(url=args.api_url, key=args.api_key, timeout=args.timeout)
    # check authentication
    client.check_authentication()
    # run program for ERA5
    ecmwf_reanalysis_retrieve(
        args.directory,
        client,
        args.year,
        SURFACE=args.surface,
        LEVEL=args.level,
        INVARIANT=args.invariant,
        MODE=args.mode,
    )
    # close connection with ECMWF datastores api server
    client = None


# run main program
if __name__ == '__main__':
    main()
