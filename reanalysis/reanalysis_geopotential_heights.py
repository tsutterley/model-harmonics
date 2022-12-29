#!/usr/bin/env python
u"""
reanalysis_geopotential_heights.py
Written by Tyler Sutterley (12/2022)
Reads temperature and specific humidity data to calculate geopotential height
    and pressure difference fields at half levels from reanalysis

INPUTS:
    Reanalysis model to run
    ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
    ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
    MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -Y X, --year X: years to run
    -V, --verbose: Output information for each output file
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use input files to define command line arguments
        added check for ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 03/2021: automatically update years to run based on current time
    Updated 01/2021: read from netCDF4 file in slices to reduce memory load
    Updated 12/2020: using argparse to set command line options
    Updated 01/2020: outputs variables as 32-bit floats instead of 64-bit floats
        clear variables and iterate years to reduce required memory
        iterate over time variable to calculate heights using incomplete files
    Updated 08/2019: adjust time scale variable for MERRA-2
    Updated 07/2018: added parameters for ERA5
    Written 03/2018
"""
from __future__ import print_function

import sys
import os
import re
import time
import logging
import netCDF4
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: reads temperature and specific humidity data to calculate
# geopotential height fields at half levels from reanalysis
def reanalysis_geopotential_heights(base_dir, MODEL, YEAR=None,
    MODE=0o775):

    # directory setup
    ddir = os.path.join(base_dir,MODEL)
    # set model specific parameters
    if (MODEL == 'ERA-Interim'):
        # invariant parameters file
        input_invariant_file = 'ERA-Interim-Invariant-Parameters.nc'
        # coordinate parameters file
        input_coordinate_file = 'ERA-Interim_coordvars.nc'
        # surface pressure file format
        input_pressure_file = 'ERA-Interim-Monthly-SP-{0:4d}.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA\-Interim\-Monthly\-Levels\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA-Interim-GPH-Levels-{0:4d}.nc'
        SURFNAME = 'z'
        ZNAME = 'z'
        VARNAME = 'sp'
        TNAME = 't'
        QNAME = 'q'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ANAME,BNAME = ('a_model_alt','b_model_alt')
        AINTERFACE,BINTERFACE = ('a_interface','b_interface')
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        UNITS = 'm**2 s**-2'
        GRAVITY = 1.0
    elif (MODEL == 'ERA5'):
        # invariant parameters file
        input_invariant_file = 'ERA5-Invariant-Parameters.nc'
        # coordinate parameters file
        input_coordinate_file = 'ERA5_coordvars.nc'
        # surface pressure file format
        input_pressure_file = 'ERA5-Monthly-SP-{0:4d}.nc'
        # regular expression pattern for finding files
        regex_pattern = r'ERA5\-Monthly\-Levels\-({0})\.nc$'
        # output file format
        output_file_format = 'ERA5-GPH-Levels-{0:4d}.nc'
        SURFNAME = 'z'
        ZNAME = 'z'
        VARNAME = 'sp'
        TNAME = 't'
        QNAME = 'q'
        DIFFNAME = 'dp'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        LEVELNAME = 'lvl'
        ANAME,BNAME = ('a_half','b_half')
        AINTERFACE,BINTERFACE = ('a_interface','b_interface')
        # hours since 1900-01-01 00:00:0.0
        TIME_LONGNAME = 'Time'
        UNITS = 'm**2 s**-2'
        GRAVITY = 1.0
    elif (MODEL == 'MERRA-2'):
        # invariant parameters file
        input_invariant_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        # coordinate parameters file
        input_coordinate_file = 'MERRA2_101.Coords_Nx.00000000.nc'
        # regular expression pattern for finding files
        regex_pattern = r'MERRA2_(\d+).instM_3d_ana_Nv.({0})(\d{{2}}).SUB.nc$'
        # output file format
        output_file_format='MERRA2_{0:0.0f}.GPH_levels.{1:4.0f}{2:02.0f}.SUB.nc'
        SURFNAME = 'PHIS'
        ZNAME = 'PHIS'
        VARNAME = 'PS'
        TNAME = 'T'
        QNAME = 'QV'
        DIFFNAME = 'dP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        LEVELNAME = 'lev'
        ANAME,BNAME = ('a_half','b_half')
        AINTERFACE,BINTERFACE = ('a_interface','b_interface')
        # minutes since start of file
        TIME_LONGNAME = 'Time'
        UNITS = 'm+2 s-2'
        GRAVITY = 1.0

    # read model orography for dimensions
    geopotential,lon,lat=ncdf_invariant(os.path.join(ddir,input_invariant_file),
        LONNAME,LATNAME,SURFNAME)
    # read parameters for calculating pressures at levels
    lev,A,B,AI,BI=ncdf_coordinates(os.path.join(ddir,input_coordinate_file),
        LEVELNAME,ANAME,BNAME,AINTERFACE,BINTERFACE)
    # Gas constant for dry air
    R_dry = 287.06

    # read each reanalysis pressure field for each year
    regex_years = r'\d{4}' if (YEAR is None) else '|'.join(map(str,YEAR))
    rx = re.compile(regex_pattern.format(regex_years), re.VERBOSE)
    input_files = [fi for fi in os.listdir(ddir) if rx.match(fi)]
    # for each reanalysis file
    for fi in sorted(input_files):
        # read input temperature and specific humidity data
        fid1 = netCDF4.Dataset(os.path.join(ddir,fi),'r')
        # extract shape from temperature variable
        ntime,nlevels,nlat,nlon = fid1.variables[TNAME].shape
        # invalid value
        fill_value = fid1.variables[TNAME]._FillValue
        # save output variables into a python dictionary.
        dinput = {}
        dinput[ZNAME] = np.zeros((ntime,nlevels,nlat,nlon),dtype=np.float32)
        dinput[DIFFNAME] = np.zeros((ntime,nlevels,nlat,nlon),dtype=np.float32)
        # model levels in reverse order
        dinput[LEVELNAME] = lev[::-1].copy()
        # extract time and time units
        dinput[TIMENAME] = np.copy(fid1.variables[TIMENAME][:])
        TIME_UNITS = fid1.variables[TIMENAME].units
        dinput[LONNAME] = lon.copy()
        dinput[LATNAME] = lat.copy()

        if MODEL in ('MERRA-2'):
            # extract date from monthly files
            MOD,YEAR,MONTH = np.array(rx.findall(fi).pop(), dtype=np.float64)
            # output monthly filename
            FILE = os.path.join(ddir,output_file_format.format(MOD,YEAR,MONTH))
            # read surface pressure
            surface_pressure = np.copy(fid1.variables[VARNAME][:])
        elif MODEL in ('ERA-Interim','ERA5'):
            # extract year from file name
            YEAR, = np.array(rx.findall(fi),dtype=np.int64)
            # output yearly filename
            FILE = os.path.join(ddir,output_file_format.format(YEAR))
            # read input surface pressure data
            pressure_file = input_pressure_file.format(YEAR)
            with netCDF4.Dataset(os.path.join(ddir,pressure_file),'r') as fid2:
                surface_pressure = np.copy(fid2.variables[VARNAME][:])

        # iterate over dates
        for t in range(ntime):
            # check dimensions for expver slice
            if (fid1.variables[VARNAME].ndim == 5):
                t_time,q_time = ncdf_expver(fid1,t,TNAME,QNAME)
            else:
                # temperature and specific humidity
                # reverse layers so bottom=0
                t_time = fid1.variables[TNAME][t,::-1,:,:]
                q_time = fid1.variables[QNAME][t,::-1,:,:]
            # calculate geopotential over model levels
            geopotential_height = np.empty((nlat,nlon),dtype=np.float32)
            # start with surface geopotential converted to units (m^2/s^2)
            geopotential_height[:,:] = geopotential*GRAVITY
            # Integrate the model layers in the atmosphere
            for k in range(nlevels):
                # calculate virtual temperature
                virtual_temp = (1.0 + 0.609133*q_time[k,:,:])*t_time[k,:,:]
                # calculate numerator and denominator for pressure ratio
                Pnum = A[k] + B[k]*surface_pressure[t,:,:]
                if ((k+1) == nlevels):
                    Pdom = 0.1
                else:
                    Pdom = A[k+1] + B[k+1]*surface_pressure[t,:,:]
                # add level to geopotential_levels
                geopotential_height[:,:] += R_dry*virtual_temp*np.log(Pnum/Pdom)
                # save level to output variable and convert to output units
                dinput[ZNAME][t,k,:,:] = geopotential_height/GRAVITY
                # calculate pressure difference between levels (at interfaces)
                Plower = AI[k] + BI[k]*surface_pressure[t,:,:]
                Pupper = AI[k+1] + BI[k+1]*surface_pressure[t,:,:]
                dinput[DIFFNAME][t,k,:,:] = Pupper - Plower

        # save to file
        ncdf_geopotential_write(dinput, fill_value, FILENAME=FILE, ZNAME=ZNAME,
            LEVELNAME=LEVELNAME, DIFFNAME=DIFFNAME, LONNAME=LONNAME,
            LATNAME=LATNAME, TIMENAME=TIMENAME, TIME_UNITS=TIME_UNITS,
            TIME_LONGNAME=TIME_LONGNAME, UNITS=UNITS)
        # set the permissions level of the output file to MODE
        os.chmod(FILE, MODE)
        # clear dinput dictionary variable
        dinput = None
        # close the input netCDF4 file
        fid1.close()

# PURPOSE: Compute the Specific Humidity from parameters (Bolton 1980)
# http://cires1.colorado.edu/~voemel/vp.html
# https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
# https://github.com/NCAR/ncl/blob/master/ni/src/lib/nfpfort/mixhum_ptrh.f
def calculate_specific_humidity(P, T, RH):
    # ratio of the molecular weights of water vapor to dry air
    epsilon = 0.622
    # calibration pressure and temperature
    pc = 6.112
    tc = 243.5
    # calculate the saturation vapor pressure in mb
    Es = pc * np.exp((17.67 * T)/(T + tc))
    # calculate the vapor pressure in mb
    Ev = Es * (RH/100.0)
    # calculate the dew point temperature
    Td = np.log(Ev/pc) * tc/(17.67 - np.log(Ev/pc))
    # calculate the specific humidity
    Q = (epsilon * Ev)/(P/100.0 - (0.378 * Ev))
    return (Q,Td)

# PURPOSE: extract temperature and specific humidity variables
# from a 5d netCDF4 dataset
# ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
def ncdf_expver(fileID, slice, TNAME, QNAME):
    ntime,nexp,nlevel,nlat,nlon = fileID.variables[TNAME].shape
    fill_value = fileID.variables[TNAME]._FillValue
    # reduced temperature and specific humidity for time
    temperature = np.ma.zeros((nlevel,nlat,nlon))
    temperature.fill_value = fill_value
    humidity = np.ma.zeros((nlevel,nlat,nlon))
    humidity.fill_value = fill_value
    # iterate over expver slices to find valid outputs
    for j in range(nexp):
        # check if any are valid for expver
        if np.any(fileID.variables[TNAME][slice,j,:,:,:]):
            # reverse layers so bottom=0
            temperature[:,:,:] = fileID.variables[TNAME][slice,j,::-1,:,:]
            humidity[:,:,:] = fileID.variables[QNAME][slice,j,::-1,:,:]
    # update mask variables
    temperature.mask = (temperature.data == temperature.fill_value)
    humidity.mask = (humidity.data == humidity.fill_value)
    # return the reduced temperature and specific humidity variables
    return (temperature,humidity)

# PURPOSE: read reanalysis invariant parameters (geopotential,lat,lon)
def ncdf_invariant(FILENAME,LONNAME,LATNAME,ZNAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        geopotential = fileID.variables[ZNAME][:].squeeze()
        longitude = fileID.variables[LONNAME][:].copy()
        latitude = fileID.variables[LATNAME][:].copy()
    return (geopotential,longitude,latitude)

# PURPOSE: read reanalysis coordinate parameters
# reverse order to go from surface to top-of-atmosphere
def ncdf_coordinates(FILENAME,LEVELNAME,ANAME,BNAME,AINTERFACE,BINTERFACE):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        levels = fileID.variables[LEVELNAME][:].copy()
        A = fileID.variables[ANAME][::-1].copy()
        B = fileID.variables[BNAME][::-1].copy()
        AI = fileID.variables[AINTERFACE][::-1].copy()
        BI = fileID.variables[BINTERFACE][::-1].copy()
    return (levels,A,B,AI,BI)

# PURPOSE: write output geopotential fields data to file
def ncdf_geopotential_write(dinput, fill_value, FILENAME=None, ZNAME=None,
    DIFFNAME=None, LEVELNAME=None, LONNAME=None, LATNAME=None, TIMENAME=None,
    TIME_UNITS=None, TIME_LONGNAME=None, UNITS=None):
    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # Defining the NetCDF dimensions
    for key in [LONNAME,LATNAME,TIMENAME,LEVELNAME]:
        fileID.createDimension(key, len(dinput[key]))

    # defining the NetCDF variables
    nc = {}
    nc[LATNAME]=fileID.createVariable(LATNAME,dinput[LATNAME].dtype,(LATNAME,))
    nc[LONNAME]=fileID.createVariable(LONNAME,dinput[LONNAME].dtype,(LONNAME,))
    nc[TIMENAME]=fileID.createVariable(TIMENAME,dinput[TIMENAME].dtype,(TIMENAME,))
    nc[LEVELNAME]=fileID.createVariable(LEVELNAME,dinput[LEVELNAME].dtype,(LEVELNAME,))
    nc[DIFFNAME] = fileID.createVariable(DIFFNAME, dinput[DIFFNAME].dtype,
        (TIMENAME,LEVELNAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)
    nc[ZNAME] = fileID.createVariable(ZNAME, dinput[ZNAME].dtype,
        (TIMENAME,LEVELNAME,LATNAME,LONNAME,), fill_value=fill_value, zlib=True)
    # filling NetCDF variables
    for key,val in dinput.items():
        nc[key][:] = dinput[key].copy()
        dinput[key] = None

    # Defining attributes for longitude and latitude
    nc[LONNAME].long_name = 'Longitude'
    nc[LONNAME].units = 'degrees_east'
    nc[LATNAME].long_name = 'Latitude'
    nc[LATNAME].units = 'degrees_north'
    # Defining attributes for time
    nc[TIMENAME].units = TIME_UNITS
    nc[TIMENAME].long_name = TIME_LONGNAME
    # Definining attributes for model levels
    nc[LEVELNAME].long_name = 'Model_Level_Number'
    # Defining attributes for geopotential height
    nc[ZNAME].long_name = 'Geopotential_Heights_on_Model_Levels'
    nc[ZNAME].units = UNITS
    # Defining attributes for pressure differences
    nc[DIFFNAME].long_name = 'Pressure_Differences_between_Levels'
    nc[DIFFNAME].units = 'Pa'

    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output NetCDF structure information
    logging.info(os.path.basename(FILENAME))
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()
    # clear nc dictionary variable
    nc = None

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads temperature and specific humidity data
            to calculate geopotential height and pressure difference
            fields at half levels from reanalysis
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = gravtk.utilities.convert_arg_line_to_args
    # command line parameters
    choices = ['ERA-Interim','ERA5','MERRA-2']
    parser.add_argument('model',
        type=str, nargs='+',
        default=['ERA5','MERRA-2'], choices=choices,
        help='Reanalysis Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to run
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Years of model outputs to run')
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
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # for each reanalysis model
    for MODEL in args.model:
        # run program
        reanalysis_geopotential_heights(args.directory, MODEL, YEAR=args.year,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
