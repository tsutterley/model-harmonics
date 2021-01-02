#!/usr/bin/env python
u"""
cds_reanalysis_retrieve.py (12/2020)
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

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: Year to retrieve
    -I, --invariant: Retrieve the model invariant parameters
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    cdsapi: Python client libraries for the CDS Web API
        https://pypi.org/project/cdsapi/

UPDATE HISTORY:
    Updated 12/2020: using argparse to set parameters
    Forked 01/2020 from ecmwf_reanalysis_retrieve.py
    Updated 07/2018: close the server connection after completion of program
        added 2-metre temperature (T2m) field as output
    Updated 03/2018: added option --mode to set the permissions level
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import cdsapi
import argparse

#-- PURPOSE: retrieve ERA5 level data for a set of years from CDS server
def cds_reanalysis_retrieve(base_dir, YEAR, INVARIANT=True, MODE=0o775):
    #-- open connection with CDS server
    server = cdsapi.Client()
    #-- parameters for ERA5 dataset
    MODEL = 'ERA5'
    model_class = "ea"
    model_dataset = "era5"
    model_grid = "0.25/0.25"
    model_levelist = "/".join(['{0:d}'.format(l) for l in range(1,137+1)])
    model_invariant_date = "2010-01-01"
    #-- output filename structure
    output_filename = "{0}-Monthly-{1}-{2:4d}.nc"
    #-- setup output directory and recursively create if currently non-existent
    ddir = os.path.join(base_dir,MODEL)
    os.makedirs(ddir, MODE) if not os.access(ddir, os.F_OK) else None

    #-- for each year
    for y in YEAR:
        #-- monthly dates to retrieve
        d="/".join(['{0:4d}{1:02d}{2:02d}'.format(y,m+1,1) for m in range(12)])

        #-- retrieve the 2-metre temperature field
        output_temperature_file = output_filename.format(MODEL,"T2m",y)
        server.retrieve("reanalysis-era5-complete", {
            "class": model_class,
            "dataset": model_dataset,
            "date": d,
            "expver": "1",
            "grid": model_grid,
            "levtype": "sfc",
            "param": "167.128",
            "stream": "moda",
            "type": "an",
            "format" : "netcdf",
        }, os.path.join(ddir,output_temperature_file))
        #-- change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_temperature_file), MODE)

        #-- retrieve the surface pressure field
        output_surface_file = output_filename.format(MODEL,"SP",y)
        server.retrieve("reanalysis-era5-complete", {
            "class": model_class,
            "dataset": model_dataset,
            "date": d,
            "expver": "1",
            "grid": model_grid,
            "levtype": "sfc",
            "param": "134.128",
            "stream": "moda",
            "type": "an",
            "format" : "netcdf",
        }, os.path.join(ddir,output_surface_file))
        #-- change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_surface_file), MODE)

        #-- retrieve model geopotential and pressure at level 1
        output_geopotential_file = output_filename.format(MODEL,"GPH",y)
        server.retrieve("reanalysis-era5-complete", {
            "class": model_class,
            "dataset": model_dataset,
            "date": d,
            "expver": "1",
            "grid": model_grid,
            "levelist": "1",
            "levtype": "ml",
            "param": "129.128/152.128",
            "stream": "moda",
            "type": "an",
            "format" : "netcdf",
        }, os.path.join(ddir,output_geopotential_file))
        #-- change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_geopotential_file), MODE)

        #-- retrieve model temperature and specific humidity
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
        #-- change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_level_file), MODE)

    #-- if retrieving the model invariant parameters
    if INVARIANT:
        output_invariant_file = '{0}-Invariant-Parameters.nc'.format(MODEL)
        server.retrieve("reanalysis-era5-complete", {
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
        }, os.path.join(ddir,output_invariant_file))
        #-- change the permissions mode to MODE
        os.chmod(os.path.join(ddir,output_invariant_file), MODE)

    #-- close connection with CDS server
    server = None

#-- Main program that calls cds_reanalysis_retrieve()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Retrieves ERA5 reanalysis netCDF4 datasets
            from the CDS Web API
            """
    )
    #-- command line parameters
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years to retrieve
    parser.add_argument('--year','-Y',
        type=int, default=range(2000,2021), nargs='+',
        help='Model years to retrieve')
    #-- retrieve the model invariant parameters
    parser.add_argument('--invariant','-I',
        default=False, action='store_true',
        help='Retrieve model invariant parameters')
    #-- permissions mode of the directories and files retrieved
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files retrieved')
    args = parser.parse_args()

    #-- run program for ERA5
    cds_reanalysis_retrieve(args.directory,args.year,
        INVARIANT=args.invariant,MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
