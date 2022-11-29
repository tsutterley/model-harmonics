#!/usr/bin/env python
u"""
cds_reanalysis_retrieve.py (11/2022)
Retrieves ERA5 reanalysis netCDF4 datasets from the CDS Web API
https://cds.climate.copernicus.eu/user/register
https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products
https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+geopotential+on+model+levels

levtype: surface (sfc), pressure levels (levtype:pl), model levels (ml)
type: an (analysis), fc (forecast)
Available parameters: http://apps.ecmwf.int/codes/grib/param-db

Get UID and API key from the CDS portal
https://cds.climate.copernicus.eu/user
Write credentials into the .cdsapirc configuration file
Alternatively write as CDSAPI_URL and CDSAPI_KEY environmental variables

COMMAND LINE OPTIONS:
    -U X, --api-url: CDS api url
    -K X, --api-key: CDS api key
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
    cdsapi: Python client libraries for the CDS Web API
        https://pypi.org/project/cdsapi/

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: added option to retrieve specific surface variables
    Updated 07/2021: added option for retrieving the model level variables
    Updated 05/2021: added option for connection timeout (in seconds)
    Updated 03/2021: added mean sea level pressure (msl) field as output
        use netCDF4 variable names for surface and invariant outputs
        automatically update years to run based on current time
    Updated 01/2021: added command line options for CDS api credentials
    Updated 12/2020: using argparse to set parameters
    Forked 01/2020 from ecmwf_reanalysis_retrieve.py
    Updated 07/2018: close the server connection after completion of program
        added 2-metre temperature (T2m) field as output
    Updated 03/2018: added option --mode to set the permissions level
    Written 03/2018
"""
from __future__ import print_function

import os
import time
import cdsapi
import argparse

# PURPOSE: retrieve ERA5 level data for a set of years from CDS server
def cds_reanalysis_retrieve(base_dir, server, YEAR,
    SURFACE=[],
    LEVEL=False,
    INVARIANT=True,
    MODE=0o775):
    # parameters for ERA5 dataset
    MODEL = 'ERA5'
    model_class = "ea"
    model_dataset = "era5"
    model_grid = "0.25/0.25"
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
        'evaporation'
        ]
    # model levels
    model_levelist = "/".join([f'{l:d}' for l in range(1,137+1)])
    # output filename structure
    output_filename = "{0}-Monthly-{1}-{2:4d}.nc"
    # setup output directory and recursively create if currently non-existent
    ddir = os.path.join(base_dir,MODEL)
    os.makedirs(ddir, MODE) if not os.access(ddir, os.F_OK) else None

    # for each year
    for y in YEAR:
        # months to retrieve
        months = [f'{m+1:02d}' for m in range(12)]
        # monthly dates to retrieve
        d = "/".join([f'{y:4d}{m}{1:02d}' for m in months])

        # for each surface variable to retrieve
        for surf in SURFACE:
            output_surface_file = output_filename.format(MODEL,surf,y)
            server.retrieve('reanalysis-era5-single-levels-monthly-means', {
                "year": str(y),
                'month': months,
                'time': '00:00',
                "grid": model_grid,
                'variable': surface_variable_dict[surf],
                "format" : "netcdf",
                'product_type': 'monthly_averaged_reanalysis',
            }, os.path.join(ddir,output_surface_file))
            # change the permissions mode to MODE
            os.chmod(os.path.join(ddir,output_surface_file), MODE)

        # if retrieving the model level data
        if LEVEL:
            # retrieve model temperature and specific humidity
            output_level_file = output_filename.format(MODEL,"Levels",y)
            server.retrieve("reanalysis-era5-complete", {
                "class": model_class,
                "dataset": model_dataset,
                "date": d,
                "expver": "1",
                "grid": model_grid,
                "levelist": model_levelist,
                "levtype": "ml",
                "param": "130.128/133.128",
                "stream": "moda",
                "type": "an",
                "format" : "netcdf",
            }, os.path.join(ddir,output_level_file))
            # change the permissions mode to MODE
            os.chmod(os.path.join(ddir,output_level_file), MODE)

    # if retrieving the model invariant parameters
    if INVARIANT:
        output_invariant_file = f'{MODEL}-Invariant-Parameters.nc'
        server.retrieve('reanalysis-era5-single-levels-monthly-means', {
            "class": model_class,
            "dataset": model_dataset,
            'year': '1979',
            'month': '01',
            'time': '00:00',
            "expver": "1",
            "grid": model_grid,
            "levtype": "sfc",
            'variable': [
                'angle_of_sub_gridscale_orography',
                'anisotropy_of_sub_gridscale_orography',
                'land_sea_mask','orography',
                'slope_of_sub_gridscale_orography',
                'standard_deviation_of_filtered_subgrid_orography',
                'standard_deviation_of_orography',
            ],
            'product_type': 'monthly_averaged_reanalysis',
            "format" : "netcdf",
        }, os.path.join(ddir,output_invariant_file))
        # change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_invariant_file), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Retrieves ERA5 reanalysis netCDF4 datasets
            from the CDS Web API
            """
    )
    # command line parameters
    # CDS api credentials
    parser.add_argument('--api-url','-U',
        type=str, default=os.environ.get('CDSAPI_URL'),
        help='CDS api url')
    parser.add_argument('--api-key','-K',
        type=str, default=os.environ.get('CDSAPI_KEY'),
        help='CDS api key')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to retrieve
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Model years to retrieve')
    # retrieve model surface variables
    # MSL: mean sea level pressure field
    # SP: surface pressure field
    # T2m: 2-metre temperature field
    # P-E: Precipitation and Evaporation fields
    choices = ['MSL','SP','T2m','P-E']
    parser.add_argument('--surface','-S',
        type=str, nargs='+', choices=choices, default=['SP'],
        help='Retrieve model surface variables')
    # retrieve the model level variables
    parser.add_argument('--level','-L',
        default=False, action='store_true',
        help='Retrieve model level variables')
    # retrieve the model invariant parameters
    parser.add_argument('--invariant','-I',
        default=False, action='store_true',
        help='Retrieve model invariant parameters')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # permissions mode of the directories and files retrieved
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files retrieved')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # open connection with CDS api server
    server = cdsapi.Client(url=args.api_url, key=args.api_key,
        timeout=args.timeout)
    # run program for ERA5
    cds_reanalysis_retrieve(args.directory, server, args.year,
        SURFACE=args.surface,
        LEVEL=args.level,
        INVARIANT=args.invariant,
        MODE=args.mode)
    # close connection with CDS api server
    server = None

# run main program
if __name__ == '__main__':
    main()
