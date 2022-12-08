#!/usr/bin/env python
u"""
era5_smb_mean.py
Written by Tyler Sutterley (12/2022)
Reads monthly ERA5 datafiles to calculate multi-annual means
    of derived surface mass balance products

ERA5 conversion table for accumulated variables
https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -m X, --mean X: Year range for mean
    -F X, --format X: input and output data format
        ascii
        netcdf
        HDF5
    -M X, --mode X: Permission mode of directories and files
    -V, --verbose: Output information for each output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 12/2021: can use variable loglevels for verbose output
    Written 10/2021
"""
from __future__ import print_function

import sys
import os
import logging
import netCDF4
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: read variables from ERA5 P-E files
def read_era5_variables(era5_flux_file):
    # python dictionary of output variables
    dinput = {}
    # read each variable of interest in ERA5 flux file
    with netCDF4.Dataset(era5_flux_file, 'r') as fileID:
        # extract geolocation variables
        dinput['latitude'] = fileID.variables['latitude'][:].copy()
        dinput['longitude'] = fileID.variables['longitude'][:].copy()
        # convert time from netCDF4 units to Julian Days
        date_string = fileID.variables['time'].units
        epoch,to_secs = gravtk.time.parse_date_string(date_string)
        dinput['time'] = gravtk.time.convert_delta_time(
            to_secs*fileID.variables['time'][:],epoch1=epoch,
            epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0) + 2400000.5
        # read each variable of interest in ERA5 flux file
        for key in ['tp','e']:
                # Getting the data from each NetCDF variable of interest
            # check dimensions for expver slice
            if (fileID.variables[key].ndim == 4):
                dinput[key] = ncdf_expver(fileID, key)
            else:
                dinput[key] = np.ma.array(fileID.variables[key][:].squeeze(),
                    fill_value=fileID.variables[key]._FillValue)
            dinput[key].mask = (dinput[key].data == dinput[key].fill_value)
    # return the output variables
    return dinput

# PURPOSE: extract variable from a 4d netCDF4 dataset
# ERA5 expver dimension (denotes mix of ERA5 and ERA5T)
def ncdf_expver(fileID, VARNAME):
    ntime,nexp,nlat,nlon = fileID.variables[VARNAME].shape
    fill_value = fileID.variables[VARNAME]._FillValue
    # reduced output
    output = np.ma.zeros((ntime,nlat,nlon))
    output.fill_value = fill_value
    for t in range(ntime):
        # iterate over expver slices to find valid outputs
        for j in range(nexp):
            # check if any are valid for expver
            if np.any(fileID.variables[VARNAME][t,j,:,:]):
                output[t,:,:] = fileID.variables[VARNAME][t,j,:,:]
    # update mask variable
    output.mask = (output.data == output.fill_value)
    # return the reduced output variable
    return output

# PURPOSE: read monthly ERA5 datasets to calculate multi-annual means
def era5_smb_mean(DIRECTORY,
    RANGE=None,
    DATAFORM=None,
    VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # sign for each product to calculate total SMB
    smb_sign = {'tp':1.0,'e':-1.0}
    # output data file format and title
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # output bad value
    fill_value = -9999.0
    # output dimensions and extents
    nlat,nlon = (721,1440)

    # attributes for output files
    attributes = {}
    attributes['varname'] = 'SMB'
    attributes['units'] = 'm w.e.'
    attributes['longname'] = 'Equivalent_Water_Thickness'
    attributes['title'] = 'ERA5 Precipitation minus Evaporation'
    attributes['source'] = ', '.join(['tp','e'])
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # mean balance flZux
    era5_mean = gravtk.spatial(nlat=nlat,nlon=nlon,
        fill_value=fill_value)
    # output data and mask
    era5_mean.data = np.zeros((nlat,nlon))
    era5_mean.mask = np.zeros((nlat,nlon),dtype=bool)
    era5_mean.time = 0.0
    count = 0.0
    # for each input file
    for Y in range(int(RANGE[0]),int(RANGE[-1])+1):
        # full path for flux file
        f1 = f'ERA5-Monthly-P-E-{Y:4d}.nc'
        era5_flux_file = os.path.join(DIRECTORY,f1)
        logging.info(era5_flux_file)
        if not os.access(era5_flux_file,os.F_OK):
            raise FileNotFoundError(f'File {f1} not in file system')
        # read netCDF4 files for variables of interest
        var = read_era5_variables(era5_flux_file)
        # days per month in year
        dpm = gravtk.time.calendar_days(Y)
        # for each month of data
        for i,t in enumerate(var['time']):
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(t,
                FORMAT='tuple')
            # spatial object for each month
            dinput = gravtk.spatial(nlat=nlat,nlon=nlon,
                fill_value=fill_value)
            dinput.lat = np.copy(var['latitude'])
            dinput.lon = np.copy(var['longitude'])
            # calculate time in year decimal
            dinput.time = gravtk.time.convert_calendar_decimal(YY,
                MM,day=DD,hour=hh,minute=mm,second=ss)
            # output data and mask
            dinput.data = np.zeros((nlat,nlon))
            dinput.mask = np.zeros((nlat,nlon),dtype=bool)
            for p in ['tp','e']:
                dinput.mask |= var[p].mask[i,:,:]
            # valid indices for all variables
            indy,indx = np.nonzero(np.logical_not(dinput.mask))
            # calculate "SMB" as precipitation minus evaporation
            # multiply by number of days to get total per month
            for p in ['tp','e']:
                dinput.data[indy,indx] += dpm[i]*var[p][i,indy,indx]*smb_sign[p]
            # update mask
            dinput.update_mask()
            # add monthly fluxes to total
            era5_mean.data += dinput.data
            era5_mean.mask |= dinput.mask
            era5_mean.time += dinput.time
            era5_mean.lat = np.copy(dinput.lat)
            era5_mean.lon = np.copy(dinput.lon)
            # add to count
            count += 1.0

    # convert from totals to means
    indy,indx = np.nonzero(np.logical_not(era5_mean.mask))
    era5_mean.data[indy,indx] /= count
    era5_mean.update_mask()
    era5_mean.time /= count

    # output ERA5 mean data file
    args = (RANGE[0], RANGE[1], suffix[DATAFORM])
    FILE = 'ERA5-Mean-P-E-{0:4d}-{1:4d}.{2}'.format(*args)
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        era5_mean.to_ascii(os.path.join(DIRECTORY,FILE),
            verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        era5_mean.to_netCDF4(os.path.join(DIRECTORY,FILE),
            verbose=VERBOSE, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        era5_mean.to_HDF5(os.path.join(DIRECTORY,FILE),
            verbose=VERBOSE, **attributes)
    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,FILE), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ERA5 datafiles to calculate
            calculate multi-annual means of derived surface mass
            balance products
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    # print information about each output file
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

    # run program with parameters
    era5_smb_mean(args.directory,
        RANGE=args.mean,
        DATAFORM=args.format,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()

