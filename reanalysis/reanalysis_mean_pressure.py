#!/usr/bin/env python
u"""
reanalysis_mean_pressure.py
Written by Tyler Sutterley (12/2022)
Calculates the mean surface pressure fields from reanalysis

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    NCEP-DOE-2: https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html
    NCEP-CFSR: https://rda.ucar.edu/datasets/ds093.1/
    JRA-55: http://jra.kishou.go.jp/JRA-55/index_en.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    --mean X: start and end year for mean
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
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for files

REFERENCES:
    JP Boy and B Chao, "Precise evaluation of atmospheric loading effects on
        Earth's time-variable gravity field", Journal of Geophysical Research:
        Solid Earth, 110(B8), (2005).
        https://doi.org/10.1029/2002JB002333

    S Swenson and J Wahr, "Estimated effects of the vertical structure of
        atmospheric mass on the time-variable geoid", Journal of Geophysical
        Research: Solid Earth, 107(B9), (2002).
        https://doi.org/10.1029/2000JB000024

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
        added check for ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: using argparse to set command line options
        using time module for operations and for extracting time units
        using spatial module for operations
    Updated 08/2019: added parameters for NCEP-CFSR, time scale for MERRA-2
    Updated 07/2018: added parameters for ERA5
    Updated 03/2018: added portions to run different reanalysis model outputs
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import copy
import logging
import netCDF4
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: read atmospheric surface pressure fields and calculates yearly mean
def reanalysis_mean_pressure(base_dir, MODEL, RANGE=None,
    VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # directory setup
    ddir = os.path.join(base_dir,MODEL)
    # set model specific parameters
    if (MODEL == 'ERA-Interim'):
        # invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA\-Interim\-Monthly\-SP\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA-Interim-Mean-SP-{0:4d}-{1:4d}.nc'
        VARNAME = 'sp'
        ZNAME = 'z'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
    elif (MODEL == 'ERA5'):
        # invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA5\-Monthly\-SP\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA5-Mean-SP-{0:4d}-{1:4d}.nc'
        VARNAME = 'sp'
        ZNAME = 'z'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
    elif (MODEL == 'MERRA-2'):
        # invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # regular expression pattern for finding files
        regex_pattern = r'MERRA2_\d{{3}}.tavgM_2d_slv_Nx.({0})(\d{{2}}).SUB.nc$'
        # output file format
        output_file_format = 'MERRA2.Mean_PS.{0:4d}-{1:4d}.nc'
        VARNAME = 'PS'
        ZNAME = 'PHIS'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'NCEP-DOE-2'):
        # invariant parameters file
        input_invariant_file = 'hgt.sfc.nc'
        # regular expression pattern for finding files
        regex_pattern = r'pres.sfc.mon.mean.({0}).nc$'
        # output file format
        output_file_format = 'pres.sfc.mean.{0:4d}-{1:4d}.nc'
        VARNAME = 'pres'
        ZNAME = 'hgt'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'NCEP-CFSR'):
        # invariant parameters file
        input_invariant_file = 'hgt.gdas.nc'
        # regular expression pattern for finding files
        regex_pattern = r'pgbh.gdas.({0}).nc$'
        # output file format
        output_file_format = 'pgbh.mean.gdas.{0:4d}-{1:4d}.nc'
        VARNAME = 'PRES_L1_Avg'
        ZNAME = 'HGT_L1_Avg'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'JRA-55'):
        # invariant parameters file
        input_invariant_file = 'll125.006_gp.2000.nc'
        # regular expression pattern for finding files
        regex_pattern = r'anl_surf125\.001_pres\.({0}).nc$'
        # output file format
        output_file_format = 'anl_surf.001_pres.mean.{0:4d}-{1:4d}.nc'
        VARNAME = 'Pressure_surface'
        ZNAME = 'GP_GDS0_SFC'
        LONNAME = 'g0_lon_1'
        LATNAME = 'g0_lat_0'
        TIMENAME = 'time'

    # read model orography for dimensions
    geopotential,lon,lat=ncdf_invariant(os.path.join(ddir,input_invariant_file),
        LONNAME,LATNAME,ZNAME)
    nlat,nlon = np.shape(geopotential)
    # read each reanalysis pressure field and calculate mean
    regex_years = r'|'.join([rf'{Y:4d}' for Y in range(RANGE[0],RANGE[1]+1)])
    rx = re.compile(regex_pattern.format(regex_years))
    input_files = [fi for fi in os.listdir(ddir) if rx.match(fi)]
    # output mean pressure field
    p_mean = gravtk.spatial()
    p_mean.lon = np.copy(lon)
    p_mean.lat = np.copy(lat)
    p_mean.time = 0.0
    p_mean.data = np.zeros((nlat,nlon))
    p_mean.mask = np.zeros((nlat,nlon),dtype=bool)
    count = 0
    # for each reanalysis file
    for fi in input_files:
        # read input data
        with netCDF4.Dataset(os.path.join(ddir,fi),'r') as fileID:
            # check dimensions for expver slice
            if (fileID.variables[VARNAME].ndim == 4):
                pressure = ncdf_expver(fileID, VARNAME)
            else:
                pressure = fileID.variables[VARNAME][:].copy()
            # use output fill value
            p_mean.fill_value = fileID.variables[VARNAME]._FillValue
            # convert time to Modified Julian Days
            delta_time=np.copy(fileID.variables[TIMENAME][:])
            date_string=fileID.variables[TIMENAME].units
            epoch,to_secs = gravtk.time.parse_date_string(date_string)
            MJD = gravtk.time.convert_delta_time(delta_time*to_secs,
                epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
        # iterate over Julian days
        for t,JD in enumerate(MJD+2400000.5):
            # add to mean pressure
            p_mean.data += pressure[t,:,:]
            p_mean.mask |= (pressure[t,:,:] == p_mean.fill_value)
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')
            # convert from calendar dates to year-decimal
            p_mean.time += gravtk.time.convert_calendar_decimal(YY,
                MM,day=DD,hour=hh,minute=mm,second=ss)
            count += 1

    # calculate mean pressure by dividing by count
    indy,indx = np.nonzero(np.logical_not(p_mean.mask))
    p_mean.data[indy,indx] /= count
    p_mean.update_mask()
    p_mean.time /= np.float64(count)

    # attributes for output files
    attributes = {}
    attributes['varname'] = copy.copy(VARNAME)
    attributes['timename'] = copy.copy(TIMENAME)
    attributes['lonname'] = copy.copy(LONNAME)
    attributes['latname'] = copy.copy(LATNAME)
    attributes['units'] = 'Pa'
    attributes['longname'] = 'surface_pressure'
    attributes['title'] = f'Surface_Pressure_from_{MODEL}_Model'
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # output to file
    FILE = output_file_format.format(RANGE[0], RANGE[1])
    # netcdf (.nc)
    p_mean.to_netCDF4(os.path.join(ddir,FILE), verbose=VERBOSE, **attributes)
    # change the permissions mode of the output file to MODE
    os.chmod(os.path.join(ddir,FILE),MODE)

# PURPOSE: extract pressure variable from a 4d netCDF4 dataset
# ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
def ncdf_expver(fileID, VARNAME):
    ntime,nexp,nlat,nlon = fileID.variables[VARNAME].shape
    fill_value = fileID.variables[VARNAME]._FillValue
    # reduced surface pressure output
    pressure = np.ma.zeros((ntime,nlat,nlon))
    pressure.fill_value = fill_value
    for t in range(ntime):
        # iterate over expver slices to find valid outputs
        for j in range(nexp):
            # check if any are valid for expver
            if np.any(fileID.variables[VARNAME][t,j,:,:]):
                pressure[t,:,:] = fileID.variables[VARNAME][t,j,:,:]
    # update mask variable
    pressure.mask = (pressure.data == pressure.fill_value)
    # return the reduced pressure variable
    return pressure

# PURPOSE: read reanalysis invariant parameters (geopotential,lat,lon)
def ncdf_invariant(FILENAME,LONNAME,LATNAME,ZNAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geopotential = fileID.variables[ZNAME][:].squeeze()
        longitude = fileID.variables[LONNAME][:].copy()
        latitude = fileID.variables[LATNAME][:].copy()
    return (geopotential,longitude,latitude)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates the mean surface pressure
            fields from reanalysis
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2','NCEP-DOE-2','NCEP-CFSR','JRA-55']
    parser.add_argument('model',
        type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[2001,2002],
        help='Start and end year range for mean')
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

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_mean_pressure(args.directory, MODEL, RANGE=args.mean,
            VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
