#!/usr/bin/env python
u"""
reanalysis_monthly_pressure.py
Written by Tyler Sutterley (07/2021)
Reads daily atmospheric pressure fields from reanalysis and outputs monthly averages

INPUTS:
    Reanalysis model to run
    NCEP-DOE-2: https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: years to run
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
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 07/2021: can use input files to define command line arguments
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 12/2020: using argparse to set command line options
        using spatial module for operations
    Updated 01/2020: changed year option to be specific years to run
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import netCDF4
import argparse
import datetime
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.spatial
import gravity_toolkit.utilities as utilities

#-- PURPOSE: read atmospheric surface pressure fields and calculate monthly mean
def reanalysis_monthly_pressure(base_dir,MODEL,YEARS,VERBOSE=False,MODE=0o775):
    #-- directory setup
    ddir = os.path.join(base_dir,MODEL)
    #-- set model specific parameters
    if (MODEL == 'NCEP-DOE-2'):
        #-- regular expression pattern for finding files
        regex_pattern = 'pres.sfc.({0:4d}).nc$'
        FILL_VALUE = 'missing_value'
        #-- output file format
        output_file_format = 'pres.sfc.mon.mean.{0:4d}.nc'
        VARNAME = 'pres'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'

    #-- for each year
    for YEAR in YEARS:
        #-- days per month (check if year is a leap year)
        dpm = [0.0] + list(gravity_toolkit.time.calendar_days(YEAR))
        cumulative_days = np.cumsum(dpm)
        #-- list of spatial data
        p_mean = []
        #-- read each reanalysis pressure field and calculate mean
        rx = re.compile(regex_pattern.format(YEAR), re.VERBOSE)
        input_files = [fi for fi in os.listdir(ddir) if rx.match(fi)]
        #-- for each input file
        for fi in input_files:
            #-- read input data
            p = gravity_toolkit.spatial().from_netCDF4(os.path.join(ddir,fi),
                varname=VARNAME, timename=TIMENAME, lonname=LONNAME,
                latname=LATNAME).transpose(axes=(1,2,0))
            p.fill_value = p.attributes['data'][FILL_VALUE]
            TIME_UNITS = p.attributes['time']['units']
            TIME_LONGNAME = p.attributes['time']['long_name']
            #-- iterate over months
            for m in range(0,12):
                #-- for each day in the month
                indices = np.arange(cumulative_days[m],cumulative_days[m+1])
                try:
                    p_mean.append(p.mean(indices=indices.astype(np.int64)))
                except:
                    break
                else:
                    p_mean[m].month = m + 1

        #-- mean pressure for each month
        p_month=gravity_toolkit.spatial().from_list(p_mean).transpose(axes=(2,0,1))
        #-- convert to python dictionary for output to netCDF4
        dinput = {}
        dinput[VARNAME] = p_month.to_masked_array()
        dinput[LATNAME] = np.copy(p_month.lat)
        dinput[LONNAME] = np.copy(p_month.lon)
        dinput[TIMENAME] = np.copy(p_month.time)

        #-- save to file
        FILE = os.path.join(ddir,output_file_format.format(YEAR))
        ncdf_pressure_write(dinput, p.fill_value, FILENAME=FILE,
            VARNAME=VARNAME, LONNAME=LONNAME, LATNAME=LATNAME, TIMENAME=TIMENAME,
            TIME_UNITS=TIME_UNITS, TIME_LONGNAME=TIME_LONGNAME, VERBOSE=VERBOSE)
        #-- set the permissions level of the output file to MODE
        os.chmod(FILE, MODE)

#-- PURPOSE: write output pressure fields data to file
def ncdf_pressure_write(dinput, fill_value, FILENAME=None, VARNAME=None,
    LONNAME=None, LATNAME=None, TIMENAME=None, TIME_UNITS=None,
    TIME_LONGNAME=None, VERBOSE=False):
    #-- opening netCDF4 file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    #-- Defining the netCDF4 dimensions
    #-- defining the netCDF4 variables
    nc = {}
    for key in [LONNAME,LATNAME,TIMENAME]:
        fileID.createDimension(key, len(dinput[key]))
        nc[key] = fileID.createVariable(key, dinput[key].dtype,(key,))
    #-- defining the main netCDF4 variable
    nc[VARNAME] = fileID.createVariable(VARNAME, dinput[VARNAME].dtype,
        (TIMENAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)
    #-- filling netCDF4 variables
    for key,val in dinput.items():
        nc[key][:] = np.copy(val)

    #-- Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'latitude'
    nc[LATNAME].units = 'degrees_north'
    #-- Defining attributes for time
    nc[TIMENAME].units = TIME_UNITS
    nc[TIMENAME].long_name = TIME_LONGNAME
    #-- Defining attributes for pressure
    nc[VARNAME].units = 'Pa'
    nc[VARNAME].long_name = 'mean_surface_pressure'

    #-- Output NetCDF structure information
    if VERBOSE:
        print(os.path.basename(FILENAME))
        print(list(fileID.variables.keys()))

    #-- Closing the NetCDF file
    fileID.close()

#-- Main program that calls reanalysis_monthly_pressure()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads daily atmospheric pressure fields
            from reanalysis and outputs monthly averages
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = utilities.convert_arg_line_to_args
    #-- command line parameters
    choices = ['NCEP-DOE-2']
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['NCEP-DOE-2'], choices=choices,
        help='Reanalysis Model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args,_ = parser.parse_known_args()

    #-- for each reanalysis model
    for MODEL in args.model:
        #-- run program
        reanalysis_monthly_pressure(args.directory, MODEL, args.year,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
