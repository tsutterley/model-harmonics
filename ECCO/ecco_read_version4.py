#!/usr/bin/env python
u"""
ecco_read_version4.py
Written by Tyler Sutterley (12/2022)

Calculates monthly ocean bottom pressure anomalies from ECCO Version 4 models
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/README
https://ecco-group.org/products-ECCO-V4r4.htm
https://ecco-group.org/user-guide-v4r4.htm

Errata document for Version 4, Revision 4 Atmospheric Pressure Forcing
https://ecco-group.org/docs/ECCO_V4r4_errata.pdf

Processes the data as described in the GRACE Tellus site
    https://grace.jpl.nasa.gov/data/get-data/ocean-bottom-pressure/
The global area average of each ocean bottom pressure map is removed

NOTES:
    Bottom Pressure Potential Anomaly (p/rhonil, m^2/s^2)
        To convert to m, divide by g (g=9.81 m/s^2)
        PHIBOT is the anomaly relative to Depth * rhonil * g
        The absolute bottom pressure in Pa is:
            Depth * rhonil * g + PHIBOT * rhonil
        rhonil = 1029 kg/m^3

INPUTS:
    ECCO version 4 models
        V4r3: Version 4, Revision 3
        V4r4: Version 4, Revision 4

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -Y X, --year X: years to run
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
    h5py: Pythonic interface to the HDF5 binary data format.
        https://www.h5py.org/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: spatial data class for reading, writing and processing data
    constants.py: calculate reference parameters for common ellipsoids

REFERENCES:
    R. J. Greatbatch, "A note on the representation of steric sea level in
        models that conserve volume rather than mass", Journal of Geophysical
        Research: Oceans, 99(C6): 12767-12771, 1994.
        https://doi.org/10.1029/94JC00847

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: lower case keyword arguments to output spatial
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 03/2021: automatically update years to run based on current time
    Updated 12/2020: use argparse to set command line parameters
        using spatial module for read/write operations
        using utilities from time module
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2018: output file with average absolute ocean bottom pressure
        run with updated 2014 GEBCO gridded bathymetry data
    Written 06/2018
"""
from __future__ import print_function

import sys
import os
import re
import logging
import datetime
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read ECCO2 V4 ocean bottom pressure data and calculate monthly
# anomalies in absolute ocean bottom pressure
def ecco_read_version4(ddir, MODEL, YEARS, RANGE=None,
    DATAFORM=None, VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # input and output subdirectories
    sd1 = f'ECCO-{MODEL}'
    sd2 = f'ECCO_{MODEL}_AveRmvd_OBP'
    # recursively create subdirectory if it doesn't exist
    if (not os.access(os.path.join(ddir,sd2), os.F_OK)):
        os.makedirs(os.path.join(ddir,sd2),MODE)
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # input variable names for each model
    if (MODEL == 'V4r3'):
        VARNAME = 'PHIBOT'
        LONNAME = 'i3'
        LATNAME = 'i2'
        TIMENAME = 'tim'
    elif (MODEL == 'V4r4'):
        VARNAME = 'PHIBOT'
        LONNAME = 'i'
        LATNAME = 'j'
        TIMENAME = 'time'
    # output dimensions and extents
    nlat,nlon = (360,720)
    extent = [-179.75,179.75,-89.75,89.75]
    # grid spacing
    dlon,dlat = (0.5,0.5)
    dphi = dlon*np.pi/180.0
    dth = dlat*np.pi/180.0
    # bad value
    fill_value = -1e+10
    # model gamma and rhonil
    gamma = 9.81
    rhonil = 1029
    # get reference parameters for WGS84 ellipsoid
    ellipsoid_params = mdlhmc.constants(ellipsoid='WGS84')
    # semimajor and semiminor axes of the ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    b_axis = ellipsoid_params.b_axis

    # read depth data from ecco_depth_version4.py
    input_depth_file = os.path.join(ddir,'DEPTH.2020.720x360.nc')
    depth = gravtk.spatial().from_netCDF4(input_depth_file,
        varname='depth', date=False)

    # read mean data from ecco_mean_version4.py
    args = (MODEL, RANGE[0], RANGE[1], suffix[DATAFORM])
    mean_file = 'ECCO_{0}_OBP_MEAN_{1:4d}-{2:4d}.{3}'.format(*args)
    # remove singleton dimensions
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        obp_mean = gravtk.spatial(
            spacing=[dlon,dlat], nlat=nlat, nlon=nlon,
            extent=extent, fill_value=fill_value).from_ascii(
            os.path.join(ddir,sd1,mean_file), date=False).squeeze()
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        obp_mean = gravtk.spatial().from_netCDF4(
            os.path.join(ddir,sd1,mean_file), date=False).squeeze()
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        obp_mean = gravtk.spatial().from_HDF5(
            os.path.join(ddir,sd1,mean_file), date=False).squeeze()

    # output average ocean bottom pressure to file
    output_average_file = f'ECCO_{MODEL}_Global_Average_OBP.txt'
    fid = open(os.path.join(ddir,sd1,output_average_file),
        mode='w', encoding='utf8')

    # compile regular expression operator for finding files for years
    regex_years = r'|'.join([rf'{y:d}' for y in YEARS])
    rx1 = re.compile(rf'PHIBOT([\.\_])({regex_years})(_\d+)?.nc$')
    # find input files
    flist = [f for f in os.listdir(os.path.join(ddir,sd1)) if rx1.match(f)]
    # read each input file
    for fi in sorted(flist):
        # Open netCDF4 datafile for reading
        PHIBOT = gravtk.spatial(fill_value=np.nan).from_netCDF4(
            os.path.join(ddir,sd1,fi),verbose=VERBOSE,
            latname=LATNAME,lonname=LONNAME,timename=TIMENAME,
            varname=VARNAME).transpose(axes=(1,2,0))
        PHIBOT.replace_invalid(fill_value)
        # time within netCDF files is days since epoch
        time_string = PHIBOT.attributes['time']['units']
        epoch1,to_secs = gravtk.time.parse_date_string(time_string)
        # read ocean bottom pressure anomalies for each month
        for m,delta_time in enumerate(to_secs*PHIBOT.time):
            # convert from ocean bottom pressure anomalies to absolute
            obp = gravtk.spatial(spacing=[dlon,dlat],nlon=nlon,
                nlat=nlat,fill_value=fill_value)
            obp.data = depth.data*rhonil*gamma + PHIBOT.data[:,:,m]*rhonil
            obp.mask = (depth.mask | PHIBOT.mask[:,:,m])
            obp.update_mask()

            # will calculate and remove the area average of the model
            obp_anomaly = gravtk.spatial(fill_value=fill_value)
            # calculate dimension variables
            obp_anomaly.lon = np.arange(extent[0],extent[1]+dlon,dlon)
            obp_anomaly.lat = np.arange(extent[2],extent[3]+dlat,dlat)
            # convert grid latitude and longitude to radians
            theta = (90.0 - obp_anomaly.lat)*np.pi/180.0
            phi = obp_anomaly.lon*np.pi/180.0
            # calculate Julian day by converting to MJD and adding offset
            JD = gravtk.time.convert_delta_time(delta_time,
                epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
                scale=1.0/86400.0) + 2400000.5
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')
            # convert from calendar dates to year-decimal
            obp_anomaly.time, = gravtk.time.convert_calendar_decimal(
                YY,MM,day=DD,hour=hh,minute=mm,second=ss)

            # global area average of each ocean bottom pressure map is removed
            # (Greatbatch correction) https://doi.org/10.1029/94JC00847
            total_area = 0.0
            total_newton = 0.0
            for k in range(0, nlat):
                # Grid point areas (ellipsoidal)
                area = np.sin(theta[k]) * np.sqrt((a_axis**2)*(b_axis**2)*
                    ((np.sin(theta[k])**2) * (np.cos(phi)**2) +
                    (np.sin(theta[k])**2)*(np.sin(phi)**2)) +
                    (a_axis**4)*(np.cos(theta[k])**2))*dphi*dth
                # calculate the grid point weight in newtons
                newtons = obp.data[k,:]*area
                # finding ocean points at each lat
                if np.count_nonzero(~obp.mask[k,:]):
                    ocean_points, = np.nonzero(~obp.mask[k,:])
                    # total area
                    total_area += np.sum(area[ocean_points])
                    # total weight in newtons
                    total_newton += np.sum(newtons[ocean_points])
            # remove global area average of each OBP map
            ratio = (total_newton/total_area)
            # Calculating Departures from the mean field
            obp_anomaly.data = (obp.data - ratio) - obp_mean.data
            obp_anomaly.mask = (obp.mask | obp_mean.mask)
            obp_anomaly.update_mask()
            # output monthly absolute bottom pressure to file
            args = (obp_anomaly.time, ratio, total_area)
            fid.write('{0:10.4f} {1:21.14e} {2:21.14e}\n'.format(*args))

            # Writing output ocean bottom pressure anomaly file
            args = (MODEL, YY, MM, suffix[DATAFORM])
            FILE = 'ECCO_{0}_AveRmvd_OBP_{1:4.0f}_{2:02.0f}.{3}'.format(*args)
            output_data(obp_anomaly, MODEL, DATAFORM=DATAFORM,
                VERBOSE=VERBOSE, FILENAME=os.path.join(ddir,sd2,FILE))
            # change the permissions mode of the output file to MODE
            os.chmod(os.path.join(ddir,sd2,FILE),MODE)

    # close output file and change the permissions to MODE
    fid.close()
    os.chmod(os.path.join(ddir,sd1,output_average_file),MODE)

# PURPOSE: wrapper function for outputting data to file
def output_data(data,MODEL,FILENAME=None,DATAFORM=None,VERBOSE=False):
    # attributes for output files
    attributes = {}
    attributes['units'] = 'Pa'
    attributes['longname'] = 'pressure_at_sea_floor'
    attributes['title'] =  f'Ocean_Bottom_Pressure_from_ECCO_{MODEL}_Model'
    attributes['reference'] = f'Output from {os.path.basename(sys.argv[0])}'
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        data.to_ascii(FILENAME,verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        data.to_netCDF4(FILENAME, verbose=VERBOSE, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        data.to_HDF5(FILENAME, verbose=VERBOSE, **attributes)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO ocean bottom pressure
            data from Version 4 models and calculates monthly
            anomalies
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['V4r3','V4r4'], choices=['V4r3','V4r4'],
        help='ECCO Version 4 Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # years to run
    now = datetime.datetime.now()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.year+1),
        help='Years of model outputs to run')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[2003,2007],
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

    # for each ECCO Version 4 model
    for MODEL in args.model:
        # run program
        ecco_read_version4(args.directory, MODEL, args.year, RANGE=args.mean,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
