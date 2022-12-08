#!/usr/bin/env python
u"""
merra_smb_cumulative.py
Written by Tyler Sutterley (12/2022)
Reads MERRA-2 datafiles to calculate monthly cumulative anomalies
    in derived surface mass balance products

From tavgM_2d_int (Vertically Integrated Diagnostics) collection:
    PRECCU (convective rain)
    PRECLS (large-scale rain)
    PRECSN (snow)
    and EVAP (evaporation)
From tavgM_2d_glc (Land Ice Surface Diagnostics) collection:
    RUNOFF (runoff over glaciated land)

INPUTS:
    SMB: Surface Mass Balance
    ACCUM: Snowfall accumulation
    PRECIP: Total Precipitation
    RAINFALL: Total Rainfall
    SUBLIM: Evaporation and Sublimation
    RUNOFF: Meltwater Runoff

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
    Updated 10/2021: using python logging for handling verbose output
        add more derived products and include sublimation and condensation
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
        sort files by month as September 2020 was reprocessed
        https://daac.gsfc.nasa.gov/information/alerts
    Updated 01/2021: use argparse to set command line parameters
        using spatial module for read/write operations
        using utilities from time module
    Updated 10/2019: changing Y/N flags to True/False
    Updated 01/2017: can output different data products
    Written 11/2016
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

# PURPOSE: read variables from MERRA-2 tavgM_2d_int and tavgM_2d_glc files
def read_merra_variables(merra_flux_file, merra_ice_surface_file):
    # python dictionary of output variables
    dinput = {}
    # read each variable of interest in MERRA-2 flux file
    with netCDF4.Dataset(merra_flux_file, 'r') as fid1:
        # extract geolocation variables
        dinput['lon'] = fid1.variables['lon'][:].copy()
        dinput['lat'] = fid1.variables['lat'][:].copy()
        # convert time from netCDF4 units to Julian Days
        date_string = fid1.variables['time'].units
        epoch,to_secs = gravtk.time.parse_date_string(date_string)
        dinput['time'] = gravtk.time.convert_delta_time(
            to_secs*fid1.variables['time'][:],epoch1=epoch,
            epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0) + 2400000.5
        # read each variable of interest in MERRA-2 flux file
        for key in ['PRECCU','PRECLS','PRECSN','EVAP']:
            # Getting the data from each NetCDF variable of interest
            dinput[key] = np.ma.array(fid1.variables[key][:].squeeze(),
                fill_value=fid1.variables[key]._FillValue)
            dinput[key].mask = (dinput[key].data == dinput[key].fill_value)
    # read each variable of interest in MERRA-2 ice surface file
    with netCDF4.Dataset(merra_ice_surface_file, 'r') as fid2:
        for key in ['RUNOFF','WESNSC']:
            # Getting the data from each NetCDF variable of interest
            dinput[key] = np.ma.array(fid2.variables[key][:].squeeze(),
                fill_value=fid2.variables[key]._FillValue)
            dinput[key].mask = (dinput[key].data == dinput[key].fill_value)
    return dinput

# PURPOSE: read monthly MERRA-2 datasets to calculate cumulative anomalies
def merra_smb_cumulative(DIRECTORY, PRODUCT, RANGE=None, DATAFORM=None,
    VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # MERRA-2 product subdirectories
    P1 = 'M2TMNXINT.5.12.4'
    P2 = 'M2TMNXGLC.5.12.4'
    # MERRA-2 output cumulative subdirectory
    cumul_sub = f'{PRODUCT}.5.12.4.CUMUL.{RANGE[0]:d}.{RANGE[1]:d}'
    # make cumulative subdirectory
    if not os.access(os.path.join(DIRECTORY,cumul_sub), os.F_OK):
        os.mkdir(os.path.join(DIRECTORY,cumul_sub), MODE)

    # regular expression operator to find datafiles (and not the xml files)
    regex_pattern = r'MERRA2_(\d+).{0}.(\d{{4}})(\d{{2}}).nc4(?!.xml)'
    # sign for each product to calculate total SMB
    smb_sign = {'PRECCU':1.0,'PRECLS':1.0,'PRECSN':1.0,'EVAP':-1.0,
        'RUNOFF':-1.0,'WESNSC':1.0}
    # titles for each output data product
    merra_products = {}
    merra_products['SMB'] = 'MERRA-2 Surface Mass Balance'
    merra_products['ACCUM'] = 'MERRA-2 Snowfall accumulation'
    merra_products['PRECIP'] = 'MERRA-2 Precipitation'
    merra_products['RAINFALL'] = 'MERRA-2 Rainfall'
    merra_products['SUBLIM'] = 'MERRA-2 Evaporation and Sublimation'
    merra_products['RUNOFF'] = 'MERRA-2 Meltwater Runoff'
    # source of each output data product
    merra_sources = {}
    merra_sources['SMB'] = ['PRECCU','PRECLS','PRECSN','EVAP','RUNOFF','WESNSC']
    merra_sources['ACCUM'] = ['PRECSN','EVAP']
    merra_sources['PRECIP'] = ['PRECCU','PRECLS','PRECSN']
    merra_sources['RAINFALL'] = ['PRECCU','PRECLS']
    merra_sources['SUBLIM'] = ['EVAP','WESNSC']
    merra_sources['RUNOFF'] = ['RUNOFF']
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')
    # output bad value
    fill_value = -9999.0
    # output dimensions and extents
    nlat,nlon = (361,576)
    extent = [-180.0,179.375,-90.0,90.0]
    # grid spacing
    dlon,dlat = (0.625,0.5)

    # attributes for output files
    attributes = {}
    attributes['varname'] = copy.copy(PRODUCT)
    attributes['units'] = 'mm w.e.'
    attributes['longname'] = 'Equivalent_Water_Thickness'
    attributes['title'] = copy.copy(merra_products[PRODUCT])
    attributes['source'] = ', '.join(merra_sources[PRODUCT])
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'

    # read mean data from merra_smb_mean.py
    args=(PRODUCT, RANGE[0], RANGE[1], suffix[DATAFORM])
    mean_file='MERRA2.tavgM_2d_{0}_mean_Nx.{1:4d}-{2:4d}.{3}'.format(*args)
    # remove singleton dimensions
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        merra_mean = gravtk.spatial(spacing=[dlon,dlat],
            nlat=nlat, nlon=nlon, extent=extent).from_ascii(
            os.path.join(DIRECTORY,mean_file), date=False).squeeze()
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        merra_mean = gravtk.spatial().from_netCDF4(
            os.path.join(DIRECTORY,mean_file),
            date=False, varname=PRODUCT).squeeze()
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        merra_mean = gravtk.spatial().from_HDF5(
            os.path.join(DIRECTORY,mean_file),
            date=False, varname=PRODUCT).squeeze()

    # find years of available data
    YEARS = sorted([d for d in os.listdir(os.path.join(DIRECTORY,P1))
        if re.match(r'\d{4}',d)])
    # check that are years are available
    CHECK = [str(Y) in YEARS for Y in range(int(YEARS[0]),int(YEARS[-1])+1)]
    if not np.all(CHECK):
        raise Exception('Not all years available on file system')
    # compile regular expression operator for flux product
    rx = re.compile(regex_pattern.format('tavgM_2d_int_Nx'), re.VERBOSE)

    # monthly cumulative anomalies
    # cumulative mass anomalies calculated by removing mean balance flux
    cumul = gravtk.spatial(nlat=nlat,nlon=nlon,fill_value=fill_value)
    cumul.lat = np.copy(merra_mean.lat)
    cumul.lon = np.copy(merra_mean.lon)
    # output data and mask
    cumul.data = np.zeros((nlat,nlon))
    cumul.mask = np.copy(merra_mean.mask)
    # for each input file
    for Y in YEARS:
        # find input files for PRODUCT
        f=[f for f in os.listdir(os.path.join(DIRECTORY,P1,Y)) if rx.match(f)]
        # sort files by month
        indices = np.argsort([rx.match(f1).group(3) for f1 in f])
        f = [f[indice] for indice in indices]
        # days per month in year
        dpm = gravtk.time.calendar_days(int(Y))
        # for each monthly file
        for M,f1 in enumerate(f):
            # extract parameters from input flux file
            MOD,Y1,M1 = rx.findall(f1).pop()
            # corresponding ice surface product file
            args = (MOD,'tavgM_2d_glc_Nx',Y1,M1)
            f2 = 'MERRA2_{0}.{1}.{2}{3}.nc4'.format(*args)
            # full path for flux and ice surface files
            merra_flux_file = os.path.join(DIRECTORY,P1,Y,f1)
            merra_ice_surface_file = os.path.join(DIRECTORY,P2,Y,f2)
            if not os.access(merra_ice_surface_file,os.F_OK):
                raise FileNotFoundError(f'File {f2} not in file system')
            # read netCDF4 files for variables of interest
            var = read_merra_variables(merra_flux_file,merra_ice_surface_file)
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(var['time'],
                FORMAT='tuple')
            # calculate the total seconds in month
            seconds = dpm[M]*24.0*60.0*60.0
            # spatial object for monthly variable
            dinput = gravtk.spatial(nlat=nlat,nlon=nlon,
                fill_value=fill_value)
            dinput.lat = np.copy(var['lat'])
            dinput.lon = np.copy(var['lon'])
            # calculate time in year decimal
            dinput.time = gravtk.time.convert_calendar_decimal(YY,
                MM,day=DD,hour=hh,minute=mm,second=ss)
            # output data and mask
            dinput.data = np.zeros((nlat,nlon))
            dinput.mask = np.zeros((nlat,nlon),dtype=bool)
            for p in ['PRECCU','PRECLS','PRECSN','EVAP','RUNOFF','WESNSC']:
                dinput.mask |= var[p].mask
            # valid indices for all variables
            indy,indx = np.nonzero(np.logical_not(dinput.mask))
            if (PRODUCT == 'SMB'):
                # calculate SMB and convert from flux to monthly
                for p in ['PRECCU','PRECLS','PRECSN','EVAP','RUNOFF','WESNSC']:
                    tmp = var[p][indy,indx]*seconds*smb_sign[p]
                    dinput.data[indy,indx] += tmp
            elif (PRODUCT == 'ACCUM'):
                # calculate accumulation and convert from flux to monthly
                for p in ['PRECSN','EVAP','WESNSC']:
                    tmp = var[p][indy,indx]*seconds*smb_sign[p]
                    dinput.data[indy,indx] += tmp
            elif (PRODUCT == 'PRECIP'):
                # calculate precipitation and convert from flux to monthly
                for p in ['PRECCU','PRECLS','PRECSN']:
                    tmp = var[p][indy,indx]*seconds
                    dinput.data[indy,indx] += tmp
            elif (PRODUCT == 'RAINFALL'):
                # calculate rainfall and convert from flux to monthly
                for p in ['PRECCU','PRECLS']:
                    tmp = var[p][indy,indx]*seconds
                    dinput.data[indy,indx] += tmp
            elif (PRODUCT == 'SUBLIM'):
                # calculate sublimation and convert from flux to monthly
                for p in ['EVAP','WESNSC']:
                    tmp = var[p][indy,indx]*seconds
                    dinput.data[indy,indx] += tmp
            elif (PRODUCT == 'RUNOFF'):
                # convert runoff from flux to monthly
                for p in ['RUNOFF']:
                    tmp = var[p][indy,indx]*seconds
                    dinput.data[indy,indx] += tmp
            # update masks
            dinput.update_mask()
            # subtract mean and add to cumulative anomalies
            cumul.data += dinput.data - merra_mean.data
            cumul.mask |= dinput.mask
            cumul.time = np.copy(dinput.time)
            cumul.update_mask()
            # copy cumulative variables to output data
            output = cumul.copy()
            # output MERRA-2 cumulative data file
            args = (MOD,PRODUCT,Y1,M1,suffix[DATAFORM])
            FILE = 'MERRA2_{0}.tavgM_2d_{1}_cumul_Nx.{2}{3}.{4}'.format(*args)
            if (DATAFORM == 'ascii'):
                # ascii (.txt)
                output.to_ascii(os.path.join(DIRECTORY,cumul_sub,FILE),
                    verbose=VERBOSE)
            elif (DATAFORM == 'netCDF4'):
                # netcdf (.nc)
                output.to_netCDF4(os.path.join(DIRECTORY,cumul_sub,FILE),
                    verbose=VERBOSE, **attributes)
            elif (DATAFORM == 'HDF5'):
                # HDF5 (.H5)
                output.to_HDF5(os.path.join(DIRECTORY,cumul_sub,FILE),
                    verbose=VERBOSE, **attributes)
            # change the permissions mode
            os.chmod(os.path.join(DIRECTORY,cumul_sub,FILE), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads MERRA-2 datafiles to calculate
            monthly cumulative anomalies in derived surface
            mass balance products
            """
    )
    # command line parameters
    choices = ['SMB','ACCUM','PRECIP','RAINFALL','SUBLIM','RUNOFF']
    parser.add_argument('product',
        type=str, nargs='+', choices=choices,
        help='MERRA-2 derived product')
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

    # run program for each input product
    for PRODUCT in args.product:
        merra_smb_cumulative(args.directory, PRODUCT, RANGE=args.mean,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()

