#!/usr/bin/env python
u"""
ecco_read_realtime.py
Written by Tyler Sutterley (12/2022)

Reads 12-hour ECCO ocean bottom pressure data from JPL
Calculates monthly anomalies on an equirectangular grid
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/KalmanFilter/
    https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Smoother/

Processes the data as described in the GRACE Tellus site
    https://grace.jpl.nasa.gov/data/get-data/ocean-bottom-pressure/
The global area average of each ocean bottom pressure map is removed

INPUTS:
    ECCO Near Real-Time models
        kf080i: Kalman filter analysis
        dr080i: RTS smoother analysis

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
    spatial.py: spatial data class for reading, writing and processing data
    time.py: utilities for calculating time operations

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
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: use argparse to set command line parameters
        using spatial module for read/write operations
        using utilities from time module
    Updated 10/2019: changing Y/N flags to True/False
    Updated 06/2019: recommending kf080i for the Kalman filtered solution
    Updated 06/2018: output file with average absolute ocean bottom pressure
    Updated 03/2018: output data in pascals
    Updated 01/2018: using getopt to set parameters
    Updated 06/2016: can use dr080g model, using __future__ print option
    Updated 02/2016: updates for new kf080g, testing different mean ranges
    Updated 06/2015: code update using regular expressions and no glob
        added main definition and DATAFORM for output ascii and HDF5 formats
    Written 02/2014
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

# PURPOSE: read ECCO ocean bottom pressure data and create monthly anomalies
# on an equirectangular grid
def ecco_read_realtime(ddir, MODEL, YEARS, RANGE=None, DATAFORM=None,
    VERBOSE=False, MODE=0o775):

    # create logger for verbosity level
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # set up regular expression for finding directories to run from YEAR
    regex_years = r'|'.join([rf'{y:d}' for y in YEARS])
    rx = re.compile(rf'{MODEL}_({regex_years})', re.VERBOSE)
    # Finding subdirectories
    input_dir = sorted([sd for sd in os.listdir(ddir) if \
        (os.path.isdir(os.path.join(ddir,sd)) & bool(rx.match(sd)))])

    # bad value
    fill_value = -1e+10
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')

    # grid parameters
    a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
    flat = 1.0/298.26# flattening of the ellipsoid
    # semiminor axis of the ellipsoid
    b_axis = (1.0 -flat)*a_axis# [m]
    # output grid spacing
    LAT_MAX = 78.5
    dlon,dlat = (1.0,1.0)
    dphi = dlon*np.pi/180.0

    # Reading Mean Ocean Bottom Pressure File
    args = (MODEL,RANGE[0],RANGE[1],suffix[DATAFORM])
    mean_file = 'ECCO_{0}_OBP_MEAN_{1:4d}-{2:4d}.{3}'.format(*args)
    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        obp_mean = gravtk.spatial(spacing=[1.0,1.0],
            nlat=158, nlon=360, extent=[0.5,359.5,-LAT_MAX,LAT_MAX],
            fill_value=fill_value).from_ascii(os.path.join(ddir,mean_file),
            date=False)
    elif (DATAFORM == 'netCDF4'):
        # netcdf (.nc)
        obp_mean = gravtk.spatial().from_netCDF4(
            os.path.join(ddir,mean_file),date=False)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        obp_mean = gravtk.spatial().from_HDF5(
            os.path.join(ddir,mean_file),date=False)

    # output subdirectory for monthly datasets
    outdir = f'ECCO_{MODEL}_AveRmvd_OBP'
    # Creating subdirectory if it doesn't exist
    if (not os.access(os.path.join(ddir,outdir), os.F_OK)):
        os.mkdir(os.path.join(ddir,outdir),MODE)

    # output average ocean bottom pressure to file
    output_average_file = f'ECCO_{MODEL}_Global_Average_OBP.txt'
    fid = open(os.path.join(ddir,outdir,output_average_file),
        mode='w', encoding='utf8')

    # for each yearly subdirectory
    for i in input_dir:
        subdir = sorted([sd for sd in os.listdir(os.path.join(ddir,i)) if
            (os.path.isdir(os.path.join(ddir,i,sd)) &
            bool(re.match(r'n10day_\d+_\d+',sd)))])
        # for each subdirectory
        for j in subdir:
            # find the input file within the subdirectory
            fi = [fi for fi in os.listdir(os.path.join(ddir,i,j)) if
                bool(re.match(r'OBP_(.*?).cdf',fi))]
            # skip subdirectory if file not found
            try:
                input_file = os.path.join(ddir,i,j,fi[0])
                os.access(input_file, os.F_OK)
            except:
                continue

            # Open ECCO CDF datafile for reading
            # change order of axes to be lat/lon/time
            obp = gravtk.spatial(fill_value=fill_value).from_netCDF4(
                input_file,verbose=VERBOSE,varname='OBP').transpose(axes=(1,2,0))
            # Getting the data from each netCDF variable
            nlat,nlon,nt = obp.shape
            # Dating scheme is hours from UNIX time (1970-01-01)
            # calculate Julian day by converting to MJD and adding offset
            time_string = obp.attributes['time']['units']
            epoch1,to_secs = gravtk.time.parse_date_string(time_string)
            JD = gravtk.time.convert_delta_time(obp.time*to_secs,
                epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
                scale=1.0/86400.0) + 2400000.5
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')

            # dlat is the difference in latitude spacing
            dlat0 = np.abs(obp.lat[0:nlat-1]-obp.lat[1:nlat])
            # used a midpoint integration method to make the hemispheres
            # symmetrical for area
            dlat1 = np.zeros((nlat))
            dlat2 = np.zeros((nlat))
            dlat1[0] = 1
            dlat2[0:-1] = np.copy(dlat0)
            dlat1[1:] = np.copy(dlat0)
            dlat2[-1] = 1
            dth = np.mean(np.array([dlat1,dlat2]),axis=0)*np.pi/180.0

            # will calculate and remove the area average of the model
            # (Greatbatch correction) https://doi.org/10.1029/94JC00847
            # Latitude spacing varies in the model
            obp_anomaly = gravtk.spatial(fill_value=fill_value)
            obp_anomaly.lon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
            obp_anomaly.lat = np.arange(-LAT_MAX,LAT_MAX+dlat,dlat)
            obp_anomaly.data = np.zeros((158,360,nt),dtype=obp.data.dtype)
            obp_anomaly.mask = np.zeros((158,360,nt),dtype=bool)
            # convert from calendar dates to year-decimal
            obp_anomaly.time = gravtk.time.convert_calendar_decimal(
                YY,MM,day=DD,hour=hh,minute=mm,second=ss)
            for t in range(0, nt):
                # the global area average of each OBP map is removed
                total_area = 0.0
                total_newton = 0.0
                for k in range(0, nlat):
                    # Grid point areas (ellipsoidal)
                    theta = (90.0 - obp.lat[k])*np.pi/180.0
                    phi = obp.lon*np.pi/180.0
                    area = np.sin(theta) * np.sqrt((a_axis**2)*(b_axis**2)*
                        ((np.sin(theta)**2) * (np.cos(phi)**2) +
                        (np.sin(theta)**2)*(np.sin(phi)**2)) +
                        (a_axis**4)*(np.cos(theta)**2))*dphi*dth[k]
                    # calculate the grid point weight in newtons
                    newtons = obp.data[k,:,t]*area
                    # finding ocean points at each lat
                    ocean_points, = np.nonzero(~obp.mask[k,:,t])
                    # total area
                    total_area += np.sum(area[ocean_points])
                    # total weight
                    total_newton += np.sum(newtons[ocean_points])
                # remove global area average of each OBP map
                ratio = (total_newton/total_area)
                obp_mean_removed = np.ma.zeros((nlat,nlon))
                obp_mean_removed.data[:,:] = obp.data[:,:,t] - ratio
                obp_mean_removed.mask = np.copy(obp.mask[:,:,t])
                # output monthly absolute bottom pressure to file
                args = (obp_anomaly.time[t], ratio, total_area)
                fid.write('{0:10.4f} {1:21.14e} {2:21.14e}\n'.format(*args))

                # interpolate to equirectangular grid
                obp_interp = np.ma.zeros((158,360))
                obp_interp.mask = np.ones((158,360),dtype=bool)
                theta = (90.0 - obp.lat)*np.pi/180.0
                th  = (90.0 - obp_anomaly.lat)*np.pi/180.0
                # for each output latitude
                for j in range(158):
                    # check if there is an exact value
                    if np.any(np.isclose(obp.lat,obp_anomaly.lat[j])):
                        k, = np.nonzero(np.isclose(obp.lat,obp_anomaly.lat[j]))
                        obp_interp.data[j,:] = obp_mean_removed.data[k,:]
                        obp_interp.mask[j,:] = obp_mean_removed.mask[k,:]
                    else:
                        # interpolate using inverse distance weights
                        # calculating the indices for the original grid
                        k, = np.nonzero((theta[0:-1] >= th[j]) &
                            (theta[1:] < th[j]))
                        # calculate distance weights
                        d1, = np.arccos(np.cos(th[j])*np.cos(theta[k]) +
                            np.sin(th[j])*np.sin(theta[k]))
                        d2, = np.arccos(np.cos(th[j])*np.cos(theta[k+1]) +
                            np.sin(th[j])*np.sin(theta[k+1]))
                        W = d1 + d2
                        # calculate interpolated value using inverse weights
                        obp_interp.data[j,:] = (obp_mean_removed.data[k,:]*d2 +
                            obp_mean_removed.data[k+1,:]*d1)/W
                        obp_interp.mask[j,:] = np.squeeze(
                            obp_mean_removed.mask[k,:] |
                            obp_mean_removed.mask[k+1,:])

                # Calculating Departures from the mean field
                obp_anomaly.data[:,:,t] = obp_interp - obp_mean.data
                obp_anomaly.mask[:,:,t] = (obp_interp.mask | obp_mean.mask)
                obp_anomaly.update_mask()

            # Calculating the monthly averages
            # data files cover the first 10 days of the next year
            ind_start_year, = np.nonzero(YY == YY[0])
            uniq_months = np.unique(MM[ind_start_year])
            for t,mm in enumerate(uniq_months):
                # Calculating the monthly anomaly
                indices, = np.nonzero(MM == mm)
                obp_monthly_anomaly = obp_anomaly.mean(indices=indices)
                obp_monthly_anomaly.update_mask()
                # output to file
                args = (MODEL,YY[0],mm,suffix[DATAFORM])
                FILE='ECCO_{0}_AveRmvd_OBP_{1:4.0f}_{2:02.0f}.{3}'.format(*args)
                output_data(obp_monthly_anomaly,MODEL,DATAFORM=DATAFORM,
                    VERBOSE=VERBOSE,FILENAME=os.path.join(ddir,outdir,FILE))
                # change the permissions mode of the output file to MODE
                os.chmod(os.path.join(ddir,outdir,FILE),MODE)

    # close output file and change the permissions to MODE
    fid.close()
    os.chmod(os.path.join(ddir,outdir,output_average_file),MODE)

# PURPOSE: wrapper function for outputting data to file
def output_data(data,MODEL,FILENAME=None,DATAFORM=None,VERBOSE=False):
    # attributes for output files
    attributes = {}
    attributes['units'] = 'Pa'
    attributes['longname'] = 'Bottom_Pressure'
    attributes['title'] = f'Ocean_Bottom_Pressure_from_ECCO-JPL_{MODEL}_Model'
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
        description="""Reads 12-hour ECCO ocean bottom pressure
            data from JPL and calculates monthly anomalies
            on an equirectangular grid
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['kf080i','dr080i'], choices=['kf080i','dr080i'],
        help='ECCO Near Real-Time Model')
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

    # for each ECCO Near Real-Time model
    for MODEL in args.model:
        # run program
        ecco_read_realtime(args.directory, MODEL, args.year, RANGE=args.mean,
            DATAFORM=args.format, VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
