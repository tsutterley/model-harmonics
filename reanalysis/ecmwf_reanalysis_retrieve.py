#!/usr/bin/env python
u"""
ecmwf_reanalysis_retrieve.py (05/2023)
Retrieves reanalysis netCDF4 datasets from the ECMWF Web API
https://software.ecmwf.int/wiki/display/CKB/How+to+download+data+via+the+ECMWF+WebAPI
https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets#AccessECMWFPublicDatasets-key

levtype: surface (sfc), pressure levels (levtype:pl), model levels (ml)
type: an (analysis), fc (forecast)
Available parameters: http://apps.ecmwf.int/codes/grib/param-db

Get UID and API key from the ECMWF portal
    Login: https://apps.ecmwf.int/auth/login/
    Key Access: https://api.ecmwf.int/v1/key/
Write credentials into the .ecmwfapirc configuration file
Alternatively write as ECMWF_API_URL, ECMWF_API_KEY and ECMWF_API_EMAIL
    environmental variables

INPUTS:
    Reanalysis dataset to retrieve (ERA-Interim or ERA5)

COMMAND LINE OPTIONS:
    -U X, --api-url: ECMWF api url
    -K X, --api-key: ECMWF api key
    -E X, --api-email: ECMWF api email
    -D X, --directory X: Working data directory
    -Y X, --year X: Year to retrieve
    -L, --level: Retrieve the model level variables
    -I, --invariant: Retrieve the model invariant parameters
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    ecmwf-api-client: Python client libraries for the ECMWF Web API
        https://pypi.org/project/ecmwf-api-client/
        https://software.ecmwf.int/wiki/display/WEBAPI/Web-API+Downloads

UPDATE HISTORY:
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
from ecmwfapi import ECMWFDataServer

# PURPOSE: retrieve ECMWF level data for a set of years
def ecmwf_reanalysis_retrieve(base_dir, server, MODEL, YEAR,
    SURFACE=[], LEVEL=False, INVARIANT=True, MODE=0o775):

    # verify input data directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # setup output directory and recursively create if currently non-existent
    ddir = base_dir.joinpath(MODEL)
    ddir.mkdir(mode=MODE, parents=True, exist_ok=True)

    # parameters for each dataset
    if (MODEL == 'ERA-Interim'):
        model_class = "ei"
        model_dataset = "interim"
        model_grid = "0.75/0.75"
        model_levelist = "/".join([f'{l:d}' for l in range(1,60+1)])
        model_invariant_date = "1989-01-01"
    elif (MODEL == 'ERA5'):
        model_class = "ea"
        model_dataset = "era5"
        model_grid = "0.25/0.25"
        model_levelist = "/".join([f'{l:d}' for l in range(1,137+1)])
        model_invariant_date = "2010-01-01"

    # surface variables
    surface_param_dict = {}
    # mean sea level pressure field
    surface_param_dict['MSL'] = "151.128"
    # surface pressure field
    surface_param_dict['SP'] = "134.128"
    # 2-metre temperature field
    surface_param_dict['T2m'] = "167.128"

    # for each year
    for y in YEAR:
        # monthly dates to retrieve
        d = "/".join([f'{y:4d}{m+1:02d}{1:02d}' for m in range(12)])

        # for each surface variable to retrieve
        for surf in SURFACE:
            output_file = ddir.joinpath(f"{MODEL}-Monthly-{surf}-{y:4d}.nc")
            server.retrieve({
                "class": model_class,
                "dataset": model_dataset,
                "date": d,
                "expver": "1",
                "grid": model_grid,
                "levtype": "sfc",
                "param": surface_param_dict[surf],
                "stream": "moda",
                "type": "an",
                "format" : "netcdf",
                "target": str(output_file),
            })
            # change the permissions mode to MODE
            output_file.chmod(mode=MODE)

        # if retrieving the model level data
        if LEVEL:
            # retrieve model temperature and specific humidity
            output_file = ddir.joinpath(f"{MODEL}-Monthly-Levels-{y:4d}.nc")
            server.retrieve({
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
                "target": str(output_file),
            })
            # change the permissions mode to MODE
            output_file.chmod(mode=MODE)

    # if retrieving the model invariant parameters
    if INVARIANT:
        output_file = f'{MODEL}-Invariant-Parameters.nc'
        server.retrieve({
            "class": model_class,
            "dataset": model_dataset,
            "date": model_invariant_date,
            "expver": "1",
            "grid": model_grid,
            "levtype": "sfc",
            "param": ("27.128/28.128/29.128/30.128/74.128/129.128/160.128/"
                "161.128/162.128/163.128/172.128"),
            "step": "0",
            "stream": "oper",
            "time": "12:00:00",
            "type": "an",
            "format" : "netcdf",
            "target": str(output_file),
        })
        # change the permissions mode to MODE
        output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Retrieves reanalysis netCDF4 datasets
            from the ECMWF Web API
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+', metavar='MODEL',
        default=['ERA-Interim'], choices=['ERA-Interim','ERA5'],
        help='Reanalysis model to retrieve')
    # ECMWF api credentials
    parser.add_argument('--api-url','-U',
        type=str, default=os.environ.get('ECMWF_API_URL'),
        help='ECMWF api url')
    parser.add_argument('--api-key','-K',
        type=str, default=os.environ.get('ECMWF_API_KEY'),
        help='ECMWF api key')
    parser.add_argument('--api-email','-E',
        type=str, default=os.environ.get('ECMWF_API_EMAIL'),
        help='ECMWF api email')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
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
    choices = ['MSL','SP','T2m']
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

    # open connection with ECMWF server
    server = ECMWFDataServer(url=args.api_url, key=args.api_key,
        email=args.api_email)
    # run program for model
    for model in args.model:
        ecmwf_reanalysis_retrieve(args.directory, server, model,
            args.year, SURFACE=args.surface, LEVEL=args.level,
            INVARIANT=args.invariant, MODE=args.mode)
    # close connection with ECMWF server
    server = None

# run main program
if __name__ == '__main__':
    main()
