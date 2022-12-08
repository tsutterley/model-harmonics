#!/usr/bin/env python
u"""
ecco_depth_version4.py
Written by Tyler Sutterley (12/2022)

Interpolates GEBCO bathymetry to ECCO Version 4 interpolated model grids
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/README
https://ecco-group.org/user-guide-v4r4.htm

GEBCO 2014/2020 Gridded bathymetry data:
https://www.bodc.ac.uk/data/hosted_data_systems/gebco_gridded_bathymetry_data/

INPUTS:
    model_file: ECCO Version 4 Model File

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
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: use argparse to set command line parameters
        using spatial module for read and write routines
    Written 10/2018
"""
from __future__ import print_function

import sys
import os
import re
import argparse
import numpy as np
import gravity_toolkit as gravtk

# PURPOSE: interpolate GEBCO bathymetry to ECCO V4 ocean model grids
def ecco_depth_version4(ddir, model_file, VERSION='2014', MODE=0o775):
    # input bathymetry model parameters
    if (VERSION == '2014'):
        FILE = os.path.join(ddir,'GEBCO_2014_2D.zip')
    elif (VERSION == '2020'):
        FILE = os.path.join(ddir,'gebco_2020_netcdf.zip')
    # read zipped file and extract file into in-memory file object
    bathymetry = gravtk.spatial().from_netCDF4(FILE,
        date=False, varname='elevation', compression='zip')
    bathymetry.data = extend_matrix(bathymetry.data,1)
    bathymetry.lon = extend_array(bathymetry.lon,1)
    bathymetry.update_mask()

    # input ECCO model parameters
    rx = re.compile(r'PHIBOT([\.\_])(\d+)(_\d+)?.nc$',re.VERBOSE)
    if rx.search(model_file).group(3):
        VARNAME,LONNAME,LATNAME,TIMENAME = ('PHIBOT','i','j','time')
    else:
        VARNAME,LONNAME,LATNAME,TIMENAME = ('PHIBOT','i3','i2','tim')
    # bad value
    fill_value = 99999.0
    # read ECCO V4 ocean model for valid points
    PHIBOT = gravtk.spatial(fill_value=np.nan).from_netCDF4(
        model_file,latname=LATNAME,lonname=LONNAME,timename=TIMENAME,
        varname=VARNAME).transpose(axes=(1,2,0)).index(0)
    PHIBOT.replace_invalid(fill_value)

    # indices of valid values
    ii,jj = np.nonzero(~PHIBOT.mask)
    # output dimensions and extents
    nlat,nlon = (360,720)
    extent = [-179.75,179.75,-89.75,89.75]
    # grid spacing
    dlon,dlat = (0.5,0.5)

    # create output data
    interp = gravtk.spatial(fill_value=fill_value)
    # calculate dimension variables
    interp.lon = np.arange(extent[0],extent[1]+dlon,dlon)
    interp.lat = np.arange(extent[2],extent[3]+dlat,dlat)
    interp.data = np.zeros((nlat,nlon))
    interp.mask = np.ones((nlat,nlon),dtype=bool)
    # iterate over indices to find valid points
    for i,j in zip(ii,jj):
        # find bathymetry points
        ilat, = np.nonzero(np.abs(interp.lat[i] - bathymetry.lat) <= 0.25)
        ilon, = np.nonzero(np.abs(interp.lon[j] - bathymetry.lon) <= 0.25)
        data_point = bathymetry.data[ilat,ilon].squeeze()
        if np.count_nonzero(data_point < 0.0):
            valid_indices, = np.nonzero(data_point <= 0.0)
            # convert from bathymetry to depth
            interp.data[i,j] = -np.mean(data_point[valid_indices])
            interp.mask[i,j] = False
    # update the mask
    interp.update_mask()

    # output file attributes
    attributes = {}
    attributes['varname'] = 'depth'
    attributes['units'] = 'm'
    attributes['longname'] = 'Depth'
    attributes['title'] = ('General Depth Chart of the Oceans, produced by the'
        'International Hydrographic Organization (IHO) and the United Nations '
        '(UNESCO) Intergovernmental Oceanographic Commission (IOC)')
    attributes['source'] = ('https://www.gebco.net/data_and_products/'
        'gridded_bathymetry_data/')
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    # output netCDF4 dataset
    bathymetry_file = f'DEPTH.{VERSION}.720x360.nc'
    interp.to_netCDF4(os.path.join(ddir,bathymetry_file),
        date=False, **attributes)
    # change the permissions mode to MODE
    os.chmod(os.path.join(ddir,bathymetry_file),MODE)

# wrapper function to extend a matrix
def extend_matrix(input_matrix,count):
    ny,nx = np.shape(input_matrix)
    temp = np.zeros((ny,nx+2*count),dtype=input_matrix.dtype)
    temp[:,0:count] = input_matrix[:,-count:]
    temp[:,count:-count] = input_matrix[:,:]
    temp[:,-count:] = input_matrix[:,0:count]
    return temp

# wrapper function to linearly extend an array
def extend_array(input_array,count):
    n = len(input_array)
    step_size = np.abs(input_array[1] - input_array[0])
    temp = np.zeros((n+2*count),dtype=input_array.dtype)
    temp[0:count] = input_array[0] - step_size*np.arange(count,0,-1)
    temp[count:-count] = input_array[:]
    temp[-count:] = input_array[-1] + step_size*np.arange(1,count+1)
    return temp

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Interpolates GEBCO bathymetry to ECCO Version 4
            interpolated model grids
            """
    )
    # command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ECCO Version 4 Model File')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # GEBCO bathymetry version year
    parser.add_argument('--version','-v',
        type=str, default='2014',
        help='GEBCO bathymetry version')
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

    # run program
    ecco_depth_version4(args.directory, args.file, VERSION=args.version,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
