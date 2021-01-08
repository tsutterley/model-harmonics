#!/usr/bin/env python
u"""
ecco_cube92_ocean_depth.py
Written by Tyler Sutterley (01/2021)

Interpolates GEBCO bathymetry to ECCO2 Cube92 ocean model grids
https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/readme.txt

GEBCO 2014/2020 Gridded bathymetry data:
https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/

INPUTS:
    ECCO2 Cube92 Model File

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -v X, --version X: GEBCO bathymetry version
    -M X, --mode X: Permission mode of directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 01/2021: use argparse to set command line parameters
        using spatial module for read and write routines
    Updated 06/2018: run with updated 2014 GEBCO gridded bathymetry data
    Written 06/2018
"""
from __future__ import print_function

import os
import re
import argparse
import numpy as np
import gravity_toolkit.spatial

#-- PURPOSE: interpolate GEBCO bathymetry to ECCO2 Cube92 ocean model grids
def ecco_cube92_ocean_depth(ddir, model_file, VERSION='2014', MODE=0o775):
    #-- input bathymetry model parameters
    if (VERSION == '2014'):
        FILE = os.path.join(ddir,'GEBCO_2014_2D.zip')
    elif (VERSION == '2020'):
        FILE = os.path.join(ddir,'gebco_2020_netcdf.zip')
    #-- read zipped file and extract file into in-memory file object
    bathymetry = gravity_toolkit.spatial().from_netCDF4(FILE,
        date=False, varname='elevation', compression='zip')
    bathymetry.data = extend_matrix(bathymetry.data,1)
    bathymetry.lon = extend_array(bathymetry.lon,1)
    bathymetry.update_mask()
    #-- bad value
    fill_value = 99999.0

    #-- read ECCO2 Cube92 ocean model for valid points
    #-- read ECCO V4 ocean model for valid points
    PHIBOT = gravity_toolkit.spatial().from_netCDF4(model_file,
        latname='LATITUDE_T',lonname='LONGITUDE_T',timename='TIME',
        varname='PHIBOT')

    #-- indices of valid values
    ii,jj = np.nonzero(~PHIBOT.mask)
    #-- output dimensions and extents
    nlat,nlon = (720,1440)
    extent = [0.125,359.875,-89.875,89.875]
    #-- grid spacing
    dlon,dlat = (0.25,0.25)

    #-- create output data
    interp = gravity_toolkit.spatial(fill_value=fill_value)
    #-- calculate dimension variables
    interp.lon = np.arange(extent[0],extent[1]+dlon,dlon)
    interp.lat = np.arange(extent[2],extent[3]+dlat,dlat)
    interp.data = np.zeros((nlat,nlon))
    interp.mask = np.ones((nlat,nlon),dtype=np.bool)
    #-- convert from 0:360 to -180:180
    gt180, = np.nonzero(PHIBOT.lon > 180.0)
    PHIBOT.lon[gt180] -= 360.0
    #-- iterate over indices to find valid points
    for i,j in zip(ii,jj):
        #-- find bathymetry points
        ilat, = np.nonzero(np.abs(PHIBOT.lat[i] - bathymetry.lat) <= 0.125)
        ilon, = np.nonzero(np.abs(PHIBOT.lon[j] - bathymetry.lon) <= 0.125)
        data_point = bathymetry.data[ilat,ilon].squeeze()
        if np.count_nonzero(data_point < 0.0):
            valid_indices, = np.nonzero(data_point <= 0.0)
            #-- convert from bathymetry to depth
            interp.data[i,j] = -np.mean(data_point[valid_indices])
            interp.mask[i,j] = False
    #-- update the mask
    interp.update_mask()

    #-- output netCDF4 dataset
    bathymetry_file = 'DEPTH.{0}.1440x720.nc'.format(VERSION)
    TITLE = ('General Depth Chart of the Oceans, produced by the'
        'International Hydrographic Organization (IHO) and the United Nations '
        '(UNESCO) Intergovernmental Oceanographic Commission (IOC)')
    REFERENCE = ('https://www.gebco.net/data_and_products/'
        'gridded_bathymetry_data/')
    interp.to_netCDF4(os.path.join(ddir,bathymetry_file), date=False,
        TITLE=TITLE, REFERENCE=REFERENCE, VARNAME='depth',
        UNITS='m', LONGNAME='Depth')
    #-- change the permissions mode to MODE
    os.chmod(os.path.join(ddir,bathymetry_file),MODE)

#-- wrapper function to extend a matrix
def extend_matrix(input_matrix,count):
    ny,nx = np.shape(input_matrix)
    temp = np.zeros((ny,nx+2*count),dtype=input_matrix.dtype)
    temp[:,0:count] = input_matrix[:,-count:]
    temp[:,count:-count] = input_matrix[:,:]
    temp[:,-count:] = input_matrix[:,0:count]
    return temp

#-- wrapper function to extend an array
def extend_array(input_array,count):
    n = len(input_array)
    step_size = np.abs(input_array[1] - input_array[0])
    temp = np.zeros((n+2*count),dtype=input_array.dtype)
    temp[0:count] = input_array[0] - step_size*np.arange(count,0,-1)
    temp[count:-count] = input_array[:]
    temp[-count:] = input_array[-1] + step_size*np.arange(1,count+1)
    return temp

#-- Main program that calls ecco_cube92_ocean_depth()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Interpolates GEBCO bathymetry to ECCO2 Cube92
            ocean model grids
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ECCO2 Cube92 Model File')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- GEBCO bathymetry version year
    parser.add_argument('--version','-v',
        type=str, default='2014',
        help='GEBCO bathymetry version')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- run program
    ecco_cube92_ocean_depth(args.directory, args.file, VERSION=args.version,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
