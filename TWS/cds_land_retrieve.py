#!/usr/bin/env python
"""
cds_land_retrieve.py (04/2025)
Retrieves ERA5-land netCDF4 datasets from the CDS Web API
https://cds.climate.copernicus.eu/user/register
https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products
https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means

Get UID and API key from the CDS portal
https://cds.climate.copernicus.eu/user
Write credentials into the .cdsapirc configuration file
Alternatively write as CDSAPI_URL and CDSAPI_KEY environmental variables

COMMAND LINE OPTIONS:
    -U X, --api-url: CDS api url
    -K X, --api-key: CDS api key
    -D X, --directory X: Working data directory
    -Y X, --year X: Year to retrieve
    -I, --invariant: Retrieve the model invariant parameters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    cdsapi: Python client libraries for the CDS Web API
        https://pypi.org/project/cdsapi/

UPDATE HISTORY:
    Written 04/2025
"""

from __future__ import print_function

import os
import time
import cdsapi
import pathlib
import logging
import argparse


# PURPOSE: retrieve ERA5-land data for a set of years from CDS client
def cds_land_retrieve(base_dir, client, YEAR, INVARIANT=True, MODE=0o775):
    # get logger
    logger = logging.getLogger(__name__)
    # verify input data directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()

    # parameters for ERA5 dataset
    MODEL = 'ERA5-Land'
    collection_id = 'reanalysis-era5-land-monthly-means'
    # setup output directory and recursively create if non-existent
    ddir = base_dir.joinpath(MODEL)
    ddir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # for each year
    for y in YEAR:
        # months to retrieve
        months = [f'{m + 1:02d}' for m in range(12)]
        output_file = ddir.joinpath(f'{MODEL}-Monthly-{y:4d}.nc')
        request = {
            'product_type': ['monthly_averaged_reanalysis'],
            'year': [str(y).zfill(4)],
            'month': months,
            'time': ['00:00'],
            'variable': [
                'snow_cover',
                'snow_depth_water_equivalent',
                'skin_reservoir_content',
                'volumetric_soil_water_layer_1',
                'volumetric_soil_water_layer_2',
                'volumetric_soil_water_layer_3',
                'volumetric_soil_water_layer_4',
            ],
            'data_format': 'netcdf',
            'download_format': 'unarchived',
        }
        # retrieve the data from the client
        client.retrieve(collection_id, request, target=str(output_file))
        # change the permissions mode to MODE
        output_file.chmod(mode=MODE)

    # if retrieving the model invariant parameters
    if INVARIANT:
        output_file = ddir.joinpath(f'{MODEL}-Invariant-Parameters.zip')
        request = {
            'product_type': ['monthly_averaged_reanalysis'],
            'year': ['1979'],
            'month': ['01'],
            'time': ['00:00'],
            'variable': [
                'glacier_mask',
                'lake_cover',
                'land_sea_mask',
                'type_of_high_vegetation',
                'type_of_low_vegetation',
            ],
            'data_format': 'netcdf',
            'download_format': 'unarchived',
        }
        # retrieve the data from the client
        client.retrieve(collection_id, request, target=str(output_file))
        # change the permissions mode to MODE
        output_file.chmod(mode=MODE)


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Retrieves ERA5-land netCDF4 datasets
            from the CDS Web API
            """
    )
    # command line parameters
    # CDS api credentials
    parser.add_argument(
        '--api-url',
        '-U',
        type=str,
        default=os.environ.get('CDSAPI_URL'),
        help='CDS api url',
    )
    parser.add_argument(
        '--api-key',
        '-K',
        type=str,
        default=os.environ.get('CDSAPI_KEY'),
        help='CDS api key',
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

    # open connection with CDS api client
    client = cdsapi.Client(
        url=args.api_url, key=args.api_key, timeout=args.timeout
    )
    # run program for ERA5-land
    cds_land_retrieve(
        args.directory,
        client,
        args.year,
        INVARIANT=args.invariant,
        MODE=args.mode,
    )
    # close connection with CDS api client
    client = None


# run main program
if __name__ == '__main__':
    main()
