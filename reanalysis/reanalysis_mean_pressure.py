#!/usr/bin/env python
u"""
reanalysis_mean_pressure.py
Written by Tyler Sutterley (02/2021)
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
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5
    time.py: utilities for calculating time operations

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
import netCDF4
import argparse
import numpy as np
import gravity_toolkit.time
import gravity_toolkit.spatial

#-- PURPOSE: read atmospheric surface pressure fields and calculates yearly mean
def reanalysis_mean_pressure(base_dir, MODEL, RANGE=None,
    VERBOSE=False, MODE=0o775):
    #-- directory setup
    ddir = os.path.join(base_dir,MODEL)
    #-- set model specific parameters
    if (MODEL == 'ERA-Interim'):
        #-- invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'ERA\-Interim\-Monthly\-SP\-({0})\.nc$'
        #-- output file format
        output_file_format = 'ERA-Interim-Mean-SP-{0:4d}-{1:4d}.nc'
        VARNAME = 'sp'
        ZNAME = 'z'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
    elif (MODEL == 'ERA5'):
        #-- invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'ERA5\-Monthly\-SP\-({0})\.nc$'
        #-- output file format
        output_file_format = 'ERA5-Mean-SP-{0:4d}-{1:4d}.nc'
        VARNAME = 'sp'
        ZNAME = 'z'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
    elif (MODEL == 'MERRA-2'):
        #-- invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        #-- regular expression pattern for finding files
        regex_pattern = 'MERRA2_\d{{3}}.tavgM_2d_slv_Nx.({0})(\d{{2}}).SUB.nc$'
        #-- output file format
        output_file_format = 'MERRA2.Mean_PS.{0:4d}-{1:4d}.nc'
        VARNAME = 'PS'
        ZNAME = 'PHIS'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'NCEP-DOE-2'):
        #-- invariant parameters file
        input_invariant_file = 'hgt.sfc.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'pres.sfc.mon.mean.({0}).nc$'
        #-- output file format
        output_file_format = 'pres.sfc.mean.{0:4d}-{1:4d}.nc'
        VARNAME = 'pres'
        ZNAME = 'hgt'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'NCEP-CFSR'):
        #-- invariant parameters file
        input_invariant_file = 'hgt.gdas.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'pgbh.gdas.({0}).nc$'
        #-- output file format
        output_file_format = 'pgbh.mean.gdas.{0:4d}-{1:4d}.nc'
        VARNAME = 'PRES_L1_Avg'
        ZNAME = 'HGT_L1_Avg'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
    elif (MODEL == 'JRA-55'):
        #-- invariant parameters file
        input_invariant_file = 'll125.006_gp.2000.nc'
        #-- regular expression pattern for finding files
        regex_pattern = r'anl_surf125\.001_pres\.({0}).nc$'
        #-- output file format
        output_file_format = 'anl_surf.001_pres.mean.{0:4d}-{1:4d}.nc'
        VARNAME = 'Pressure_surface'
        ZNAME = 'GP_GDS0_SFC'
        LONNAME = 'g0_lon_1'
        LATNAME = 'g0_lat_0'
        TIMENAME = 'time'

    #-- read model orography for dimensions
    geopotential,lon,lat=ncdf_invariant(os.path.join(ddir,input_invariant_file),
        LONNAME,LATNAME,ZNAME)
    nlat,nlon = np.shape(geopotential)
    #-- read each reanalysis pressure field and calculate mean
    regex_years = '|'.join(['{0:4d}'.format(Y) for Y in range(RANGE[0],RANGE[1]+1)])
    rx = re.compile(regex_pattern.format(regex_years))
    input_files = [fi for fi in os.listdir(ddir) if rx.match(fi)]
    #-- output mean pressure field
    p_mean = gravity_toolkit.spatial()
    p_mean.lon = np.copy(lon)
    p_mean.lat = np.copy(lat)
    p_mean.time = 0.0
    p_mean.data = np.zeros((nlat,nlon))
    p_mean.mask = np.zeros((nlat,nlon),dtype=bool)
    count = 0
    #-- for each reanalysis file
    for fi in input_files:
        #-- read input data
        with netCDF4.Dataset(os.path.join(ddir,fi),'r') as fileID:
            pressure = fileID.variables[VARNAME][:].copy()
            p_mean.fill_value = fileID.variables[VARNAME]._FillValue
            #-- convert time to Modified Julian Days
            delta_time=np.copy(fileID.variables[TIMENAME][:])
            date_string=fileID.variables[TIMENAME].units
            epoch,to_secs=gravity_toolkit.time.parse_date_string(date_string)
            MJD=gravity_toolkit.time.convert_delta_time(delta_time*to_secs,
                epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
        #-- iterate over Julian days
        for t,JD in enumerate(MJD+2400000.5):
            #-- add to mean pressure
            p_mean.data += pressure[t,:,:]
            p_mean.mask |= (pressure[t,:,:]==p_mean.fill_value)
            #-- convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(JD,
                FORMAT='tuple')
            #-- convert from calendar dates to year-decimal
            p_mean.time += gravity_toolkit.time.convert_calendar_decimal(YY,
                MM,day=DD,hour=hh,minute=mm,second=ss)
            count += 1

    #-- calculate mean pressure by dividing by count
    indy,indx = np.nonzero(~p_mean.mask)
    p_mean.data[indy,indx] /= count
    p_mean.update_mask()
    p_mean.time /= np.float(count)

    #-- output to file
    FILE = output_file_format.format(RANGE[0], RANGE[1])
    TITLE = 'Mean_Surface_Pressure_from_{0}_Model'.format(MODEL)
    #-- netcdf (.nc)
    p_mean.to_netCDF4(os.path.join(ddir,FILE), verbose=VERBOSE,
        varname=VARNAME, timename=TIMENAME, lonname=LONNAME, latname=LATNAME,
        UNITS='Pa', LONGNAME='surface_pressure', TITLE=TITLE)
    #-- change the permissions mode of the output file to MODE
    os.chmod(os.path.join(ddir,FILE),MODE)

#-- PURPOSE: read reanalysis invariant parameters (geopotential,lat,lon)
def ncdf_invariant(FILENAME,LONNAME,LATNAME,ZNAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geopotential = fileID.variables[ZNAME][:].squeeze()
        longitude = fileID.variables[LONNAME][:].copy()
        latitude = fileID.variables[LATNAME][:].copy()
    return (geopotential,longitude,latitude)

#-- Main program that calls reanalysis_mean_pressure()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Calculates the mean surface pressure
            fields from reanalysis
            """
    )
    #-- command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2','NCEP-DOE-2','NCEP-CFSR','JRA-55']
    parser.add_argument('model',
        metavar='MODEL', type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- start and end years to run for mean
    parser.add_argument('--mean',
        metavar=('START','END'), type=int, nargs=2,
        default=[2001,2002],
        help='Start and end year range for mean')
    #-- print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- for each reanalysis model
    for MODEL in args.model:
        #-- run program
        reanalysis_mean_pressure(args.directory, MODEL, RANGE=args.mean,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
