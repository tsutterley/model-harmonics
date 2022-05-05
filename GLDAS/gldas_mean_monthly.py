#!/usr/bin/env python
u"""
gldas_mean_monthly.py
Written by Tyler Sutterley (05/2022)

Reads GLDAS monthly datafiles from http://ldas.gsfc.nasa.gov/gldas/
Adding Soil Moisture, snow water equivalent (SWE) and total canopy storage
Calculates the time-averaged grid from a determined range of years

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

INPUTS:
    GLDAS land surface model
        CLM: Common Land Model (CLM)
        CLSM: Catchment Land Surface Model (CLSM)
        MOS: Mosaic model
        NOAH: Noah model
        VIC: Variable Infiltration Capacity (VIC) model

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -m X, --mean X: Year range for mean
    -S X, --spacing X: spatial resolution of models to run
        10: 1.0 degrees latitude/longitude
        025: 0.25 degrees latitude/longitude
    -v X, --version X: GLDAS model version
    -F X, --format X: Input and output data format
        ascii
        netCDF4
        HDF5
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    pygrib: Python interface for reading and writing GRIB data
        https://pypi.python.org/pypi/pygrib
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: add try/except for pygrib import
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: set spatial variables for both 025 and 10 cases
        using utilities from time module
    Updated 10/2020: use argparse to set command line parameters
    Updated 08/2020: flake8 compatible regular expression strings
    Updated 06/2020: using spatial data class for input and output operations
    Updated 02/2020: modified netCDF4 read function for CLSM and VIC version 2.1
    Updated 10/2019: can read subsetted netCDF4 files for VIC models
        changing Y/N flags to True/False
    Updated 06/2018: using python3 compatible octal and input
    Updated 03/2018: include estimates of total canopy water storage
    Forked 03/2018 from gldas_mean_3H.py
    Updated 03/2018: using getopt to set parameters.  can read netCDF4 files.
    Updated 06/2016: updated to use __future__ print function
        can output different GLDAS mean fields
    Updated 05/2015-06/2015: code update using regular expressions and no glob
        added DATAFORM for output ascii and HDF5 formats
    Updated 02/2014 with quick updates but should be fully updated
    Written 04/2013
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import argparse
import warnings
import numpy as np
from gravity_toolkit.time import convert_calendar_decimal
from gravity_toolkit.spatial import spatial
try:
    import pygrib
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("pygrib not available")

#-- GLDAS models
gldas_products = {}
gldas_products['CLM'] = 'GLDAS Common Land Model (CLM)'
gldas_products['CLSM'] = 'GLDAS Catchment Land Surface Model (CLSM)'
gldas_products['MOS'] = 'GLDAS Mosaic model'
gldas_products['NOAH'] = 'GLDAS Noah model'
gldas_products['VIC'] = 'GLDAS Variable Infiltration Capacity (VIC) model'

def gldas_mean_monthly(base_dir, MODEL, RANGE=None, SPATIAL=None, VERSION=None,
    DATAFORM=None, VERBOSE=False, MODE=0o775):

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    #-- Version flags
    V1,V2 = ('_V1','') if (VERSION == '1') else ('','.{0}'.format(VERSION))
    #-- dimensions of spatial fields from SPATIAL variable
    if (SPATIAL == '025'):
        nlon,nlat = (1440,600)
    elif (SPATIAL == '10'):
        nlon,nlat = (360,150)
    #-- subdirectory for model monthly products at spacing for version
    subdir = "GLDAS_{0}{1}_{2}{3}".format(MODEL,SPATIAL,'M',V2)
    #-- directory for GLDAS model
    ddir = os.path.join(base_dir, subdir)
    #-- years to run for mean
    year_dir = ['{0:4d}'.format(sd) for sd in range(RANGE[0],RANGE[1]+1)]
    #-- output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- compile regular expression pattern for finding files
    GLDAS_SUFFIX = r'nc4|grb|grb\.SUB\.nc4'
    regex_pattern = r'GLDAS_{0}{1}_{2}\.A(\d{{4}})(\d{{2}})\.(\d+)\.({3})$'
    rx = re.compile(regex_pattern.format(MODEL,SPATIAL,'M',GLDAS_SUFFIX))
    #-- ndates is the number of monthly measurements between date range
    ndates = 12*(RANGE[1]-RANGE[0]+1)

    #-- allocate for mean TWC and date
    twc = spatial()
    twc.data = np.zeros((nlat,nlon,ndates))
    twc.mask = np.ones((nlat,nlon,ndates),dtype=bool)
    twc.time = np.zeros((ndates))
    #-- create counter for dates
    c = 0
    #-- for each year within years_range
    for i,yr in enumerate(year_dir):
        #-- find all GRIB/netCDF4 files within directory
        f = [f for f in os.listdir(os.path.join(ddir,yr)) if rx.match(f)]
        #-- for each GRIB/netCDF4 file
        for fi in sorted(f):
            #-- Getting date information from file
            YY,MM,VF,SFX = rx.findall(fi).pop()
            #-- read GRIB or netCDF4 file
            if (SFX == 'grb'):
                SM,SWE,CW,twc.lat,twc.lon,twc.fill_value = \
                    grib_twc_read(os.path.join(ddir,yr,fi))
            elif (SFX == 'nc4'):
                SM,SWE,CW,twc.lat,twc.lon,twc.fill_value = \
                    ncdf_twc_read(os.path.join(ddir,yr,fi))
            elif (SFX == 'grb.SUB.nc4'):
                SM,SWE,CW,twc.lat,twc.lon,twc.fill_value = \
                    subset_twc_read(os.path.join(ddir,yr,fi))
            #-- converting from kg/m^2 to cm water equivalent (cmwe)
            ii,jj = np.nonzero(SWE != twc.fill_value)
            twc.data[ii,jj,c] = 0.1*(SM[ii,jj] + SWE[ii,jj] + CW[ii,jj])
            #-- set the mask for valid points
            twc.mask[ii,jj,c] = False
            #-- calculate date
            twc.time[c] = convert_calendar_decimal(np.int64(YY),np.int64(MM))
            #-- add 1 to counter
            c += 1

    #-- calculate mean TWC and date
    twc_mean = twc.mean()

    #-- output to file
    args = (MODEL,SPATIAL,RANGE[0],RANGE[1],suffix[DATAFORM])
    FILE = 'GLDAS_{0}{1}_TWC_MEAN_{2:4d}-{3:4d}.{4}'.format(*args)
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        twc_mean.to_ascii(os.path.join(ddir,FILE),date=False,
            verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netCDF4 (.nc)
        twc_mean.to_netCDF4(os.path.join(ddir,FILE),verbose=VERBOSE,
            units='cmwe',longname='Equivalent Water Thickness',
            title=gldas_products[MODEL],date=False)
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        twc_mean.to_HDF5(os.path.join(ddir,FILE),verbose=VERBOSE,
            units='cmwe',longname='Equivalent Water Thickness',
            title=gldas_products[MODEL],date=False)
    #-- change the permissions mode
    os.chmod(os.path.join(ddir,FILE), MODE)

#-- PURPOSE: read a GLDAS GRIB file for snow_water_eq and soil_moisture
def grib_twc_read(FILENAME):
    #-- Opening GRB file of year and month
    fileID = pygrib.open(FILENAME)
    #-- bad value
    fill_value = -9999.0
    #-- Getting variables
    #-- Soil Moisture:
    soil_values = []
    #-- for each soil moisture level
    nsoil = len(fileID.select(name='Soil Moisture'))
    for level in range(nsoil):
        soilinp = fileID.select(name='Soil Moisture')[level]
        #-- converting masked arrays to standard arrays masked with fill_value
        soil_level = np.ma.filled(soilinp['values'],fill_value=fill_value)
        soil_values.append(soil_level)
    #-- taking the sum of the soil moisture levels
    soil_moisture = np.sum(soil_values,axis=0)
    #-- Snow Water Equivalent:
    snowinp = fileID.select(name='Snow Fall water equivalent')[0]
    #-- converting masked arrays to standard arrays masked with fill_value
    snow_water_eq = np.ma.filled(snowinp['values'],fill_value=fill_value)
    #-- Canopy Water Amount: (only listed by indicator)
    #-- data provided by GRACE Tellus neglected Total canopy water storage
    canopyinp = fileID.select(indicatorOfParameter=71)[0]
    #-- converting masked arrays to standard arrays masked with fill_value
    canopy_water = np.ma.filled(canopyinp['values'],fill_value=fill_value)
    #-- Getting lat and lon from snow product
    lat, lon = snowinp.latlons()
    #-- lon and lat are matrices.. getting the unique arrays
    glat = lat[:,0]
    glon = lon[0,:]
    #-- return values
    return (soil_moisture,snow_water_eq,canopy_water,glat,glon,fill_value)

#-- PURPOSE: read a GLDAS netCDF4 file for snow_water_eq and soil_moisture
def ncdf_twc_read(FILENAME):
    #-- Opening netCDF4 file of year and month
    fileID = netCDF4.Dataset(FILENAME,'r')
    #-- read variables from netCDF4 file
    lon = fileID.variables['lon'][:].squeeze()
    lat = fileID.variables['lat'][:].squeeze()
    #-- read snow water equivalent (SWE)
    fill_value = fileID.variables['SWE_inst']._FillValue
    snow_water_eq = fileID.variables['SWE_inst'][:].squeeze()
    #-- read plant canopy water amount
    canopy_water = fileID.variables['CanopInt_inst'][:].squeeze()
    #-- allocate for soil moisture
    soil_moisture = np.full((len(lat),len(lon)),fill_value)
    #-- regular expression pattern for soil moisture levels
    regex = re.compile(r'SoilMoi(.*?)_inst',re.VERBOSE)
    keys = [k for k in fileID.variables.keys() if regex.match(k)]
    for key in keys:
        val = fileID.variables[key][:].squeeze()
        ii,jj = np.nonzero(val != fill_value)
        soil_moisture[ii,jj] += val[ii,jj]
    #-- close the netCDF4 file
    fileID.close()
    #-- return values
    return (soil_moisture,snow_water_eq,canopy_water,lat,lon,fill_value)

#-- PURPOSE: read subset GLDAS netCDF4 file for snow_water_eq and soil_moisture
def subset_twc_read(FILENAME):
    #-- Opening netCDF4 file of year and month
    fileID = netCDF4.Dataset(FILENAME,'r')
    #-- read variables from netCDF4 file
    lon = fileID.variables['lon'][:].squeeze()
    lat = fileID.variables['lat'][:].squeeze()
    #-- read snow water equivalent (SWE)
    fill_value = fileID.variables['SWE']._FillValue
    snow_water_eq = fileID.variables['SWE'][:].squeeze()
    #-- read plant canopy water amount
    canopy_water = fileID.variables['Canint'][:].squeeze()
    #-- sum soil moisture over model layers
    soil_moisture = np.sum(fileID.variables['SoilMoist'][:].squeeze(),axis=0)
    #-- close the netCDF4 file
    fileID.close()
    #-- return values
    return (soil_moisture,snow_water_eq,canopy_water,lat,lon,fill_value)

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads GLDAS monthly datafiles to calculate the
            multi-annual mean total water storage from soil moisture,
            snow water equivalent and total canopy storage
            """
    )
    #-- command line parameters
    parser.add_argument('model',
        type=str, nargs='+', choices=gldas_products.keys(),
        help='GLDAS land surface model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[2003,2007],
        help='Start and end year range for mean')
    #-- GLDAS model version
    parser.add_argument('--version','-v',
        type=str, default='2.1',
        help='GLDAS model version')
    #-- model spatial resolution
    #-- 10: 1.0 degrees latitude/longitude
    #-- 025: 0.25 degrees latitude/longitude
    parser.add_argument('--spacing','-S',
        type=str, default='10', choices=['10','025'],
        help='Spatial resolution of models to run')
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- for each GLDAS model
    for MODEL in args.model:
        #-- run program
        gldas_mean_monthly(args.directory, MODEL, VERSION=args.version,
            SPATIAL=args.spacing, RANGE=args.mean, DATAFORM=args.format,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
