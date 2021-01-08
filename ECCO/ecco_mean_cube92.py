#!/usr/bin/env python
u"""
ecco_mean_cube92.py
Written by Tyler Sutterley (01/2021)

Calculates mean of ocean bottom pressure data from the ECCO2 Cube92 ocean model
https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/readme.txt

Computes the time-average ocean bottom pressure map between two dates
Processes the data as described in the GRACE Tellus site
    https://grace.jpl.nasa.gov/data/get-data/ocean-bottom-pressure/
The global area average of each ocean bottom pressure map is removed

NOTES:
    Bottom Pressure Potential Anomaly (p/rhonil, m^2/s^2)
        To convert to m, divide by g (g=9.81 m/s^2)
        PHIBOT is the anomaly relative to Depth * rhonil * g
        The absolute bottom pressure in Pa is:
            Depth * rhonil * g + PHIBOT * rhonil
        rhonil = 1027.5 kg/m^3

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
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    netCDF4: Python interface to the netCDF C library
        (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
    h5py: Pythonic interface to the HDF5 binary data format.
        (https://www.h5py.org/)

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: spatial data class for reading, writing and processing data
        ncdf_read.py: reads input spatial data from netCDF4 files
        hdf5_read.py: reads input spatial data from HDF5 files
        ncdf_write.py: writes output spatial data to netCDF4
        hdf5_write.py: writes output spatial data to HDF5
    ref_ellipsoid.py: calculate reference parameters for common ellipsoids

REFERENCES:
    R. J. Greatbatch, "A note on the representation of steric sea level in
        models that conserve volume rather than mass", Journal of Geophysical
        Research: Oceans, 99(C6): 12767-12771, 1994.
        https://doi.org/10.1029/94JC00847

UPDATE HISTORY:
    Updated 01/2021: use argparse to set command line parameters
        using spatial module for read/write operations
        using utilities from time module
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2018: run with updated 2014 GEBCO gridded bathymetry data
    Written 06/2018
"""
from __future__ import print_function

import os
import re
import argparse
import numpy as np
import gravity_toolkit.spatial
import gravity_toolkit.time
from geoid_toolkit.ref_ellipsoid import ref_ellipsoid

#-- PURPOSE: read ECCO2 Cube92 ocean bottom pressure data and calculate mean
def ecco_mean_cube92(ddir, RANGE=None, DATAFORM=None, VERBOSE=False, MODE=0o775):
    #-- input and output subdirectories
    sd = 'cube92_latlon_quart_90S90N'
    #-- output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    #-- input variable names
    VARNAME = 'PHIBOT'
    LONNAME = 'LONGITUDE_T'
    LATNAME = 'LATITUDE_T'
    TIMENAME = 'TIME'
    #-- output dimensions and extents
    nlat,nlon = (720,1440)
    extent = [0.125,359.875,-89.875,89.875]
    #-- grid spacing
    dlon,dlat = (0.25,0.25)
    dphi = dlon*np.pi/180.0
    dth = dlat*np.pi/180.0
    #-- gamma and rhonil
    gamma = 9.81
    rhonil = 1027.5
    #-- get reference parameters for WGS84 ellipsoid
    WGS84 = ref_ellipsoid('WGS84')
    #-- semimajor and semiminor axes of the ellipsoid [m]
    a_axis = WGS84['a']
    b_axis = WGS84['b']

    #-- read depth data from ecco_cube92_ocean_depth.py
    input_depth_file = os.path.join(ddir,'DEPTH.2020.720x360.nc')
    depth = gravity_toolkit.spatial().from_netCDF4(input_depth_file,
        varname='depth', date=False)

    #-- compile regular expression operator for years
    year_regex = '|'.join('{0:d}'.format(y) for y in range(RANGE[0],RANGE[1]+1))
    rx1 = re.compile(r'PHIBOT\.(\d+)x(\d+)\.({0})(\d{{2}}).nc$'.format(year_regex))
    #-- find input files
    input_files = [fi for fi in os.listdir(ddir) if rx1.match(fi)]

    #-- output multi-annual mean
    obp_mean = gravity_toolkit.spatial()
    obp_mean.data = np.zeros((nlat,nlon),dtype=np.float)
    obp_mean.mask = np.zeros((nlat,nlon),dtype=np.bool)
    obp_mean.time = 0.0
    #-- calculate dimension variables
    obp_mean.lon = np.arange(extent[0],extent[1]+dlon,dlon)
    obp_mean.lat = np.arange(extent[2],extent[3]+dlat,dlat)
    #-- convert grid latitude and longitude to radians
    theta = (90.0 - obp_mean.lat)*np.pi/180.0
    phi = obp_mean.lon*np.pi/180.0
    #-- counter variable for dates
    count = 0.0
    #-- read each input file
    for t,fi in enumerate(input_files):
        #-- Open netCDF4 datafile for reading
        PHIBOT = gravity_toolkit.spatial().from_netCDF4(
            os.path.join(ddir,sd,fi),verbose=VERBOSE,
            latname=LATNAME,lonname=LONNAME,timename=TIMENAME,
            varname=VARNAME).transpose(axes=(1,2,0))
        #-- time within netCDF files is days since 1992-01-01
        time_string = PHIBOT.attributes['time']['units']
        epoch1,to_secs = gravity_toolkit.time.parse_date_string(time_string)

        #-- calculate Julian day by converting to MJD and adding offset
        JD = gravity_toolkit.time.convert_delta_time(to_secs*PHIBOT.time,
            epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
            scale=1.0/86400.0) + 2400000.5
        #-- convert from Julian days to calendar dates
        YY,MM,DD,hh,mm,ss = gravity_toolkit.time.convert_julian(JD,
            FORMAT='tuple')
        #-- convert from calendar dates to year-decimal
        obp_mean.time += gravity_toolkit.time.convert_calendar_decimal(YY,MM,
            day=DD,hour=hh,minute=mm,second=ss)

        #-- convert from ocean bottom pressure anomalies to absolute
        obp = gravity_toolkit.spatial(spacing=[dlon,dlat],nlon=nlon,
            nlat=nlat,fill_value=PHIBOT.fill_value)
        obp.data = depth.data*rhonil*gamma + PHIBOT.data[:,:,0]*rhonil
        obp.mask = (depth.mask | PHIBOT.mask[:,:,0])
        obp.update_mask()

        #-- global area average of each ocean bottom pressure map is removed
        #-- (Greatbatch correction) https://doi.org/10.1029/94JC00847
        total_area = 0.0
        total_newton = 0.0
        for k in range(0, nlat):
            #-- Grid point areas (ellipsoidal)
            area = np.sin(theta[k]) * np.sqrt((a_axis**2)*(b_axis**2)*
                ((np.sin(theta[k])**2) * (np.cos(phi)**2) +
                (np.sin(theta[k])**2)*(np.sin(phi)**2)) +
                (a_axis**4)*(np.cos(theta[k])**2))*dphi*dth
            #-- calculate the grid point weight in newtons
            newtons = obp.data[k,:]*area
            #-- finding ocean points at each lat
            if np.count_nonzero(~obp.mask[k,:]):
                ocean_points, = np.nonzero(~obp.mask[k,:])
                #-- total area
                total_area += np.sum(area[ocean_points])
                #-- total weight in newtons
                total_newton += np.sum(newtons[ocean_points])
        #-- remove global area average of each OBP map
        ratio = (total_newton/total_area)
        obp_mean.data += (obp.data - ratio)
        obp_mean.mask |= np.copy(obp.mask)
        count += 1.0

    #-- convert from totals to means
    indy,indx = np.nonzero(~obp_mean.mask)
    obp_mean.data[indy,indx] /= count
    obp_mean.update_mask()
    obp_mean.replace_invalid(PHIBOT.fill_value)
    obp_mean.time /= count

    #-- output to file
    args = (RANGE[0], RANGE[1], suffix[DATAFORM])
    FILE = 'ECCO_Cube92_OBP_MEAN_{0:4d}-{1:4d}.{2}'.format(*args)
    TITLE = 'Mean_Ocean_Bottom_Pressure_from_ECCO2_Cube92_Model'
    if (DATAFORM == 'ascii'):
        #-- ascii (.txt)
        obp_mean.to_ascii(os.path.join(ddir,sd,FILE),verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        #-- netcdf (.nc)
        obp_mean.to_netCDF4(os.path.join(ddir,sd,FILE), verbose=VERBOSE,
            UNITS='Pa', LONGNAME='pressure_at_sea_floor', TITLE=TITLE)
    elif (DATAFORM == 'HDF5'):
        #-- HDF5 (.H5)
        obp_mean.to_HDF5(os.path.join(ddir,sd,FILE), verbose=VERBOSE,
            UNITS='Pa', LONGNAME='pressure_at_sea_floor', TITLE=TITLE)
    #-- change the permissions mode of the output file to MODE
    os.chmod(os.path.join(ddir,sd,FILE),MODE)

#-- Main program that calls ecco_mean_cube92()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO2 Cube92 ocean bottom
            pressure data and calculates multi-annual means
            """
    )
    #-- command line parameters
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
    #-- input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    args = parser.parse_args()

    #-- run program
    ecco_mean_cube92(args.directory, RANGE=args.mean, DATAFORM=args.format,
        VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
