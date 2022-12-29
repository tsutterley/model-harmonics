#!/usr/bin/env python
u"""
ecco_mean_llc_tiles.py
Written by Tyler Sutterley (12/2022)

Calculates mean of tiled ocean bottom pressure data from the ECCO ocean model
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/nctiles_monthly
https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/nctiles_monthly

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
        rhonil = 1029 kg/m^3

INPUTS:
    ECCO LLC tile models
        V4r4: Version 4, Revision 4
        V5alpha: Version 5, Alpha release

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -m X, --mean X: Year range for mean
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

PROGRAM DEPENDENCIES:
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
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import re
import logging
import netCDF4
import datetime
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read ECCO tiled ocean bottom pressure data and calculate mean
def ecco_mean_llc_tiles(ddir, MODEL, RANGE=None, MODE=0o775):

    # input and output subdirectories
    DIRECTORY = os.path.join(ddir,f'ECCO-{MODEL}','nctiles_monthly')

    # input variable names for each model
    if (MODEL == 'V4r4'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        TIMENAME = 'time'
        AREANAME = 'rA'
        MASKNAME = 'maskC'
        Nt,Nj,Ni = (13,90,90)
    elif (MODEL == 'V5alpha'):
        LONNAME = 'XC'
        LATNAME = 'YC'
        ZNAME = 'Depth'
        VARNAME = 'PHIBOT'
        TIMENAME = 'time'
        AREANAME = 'rA'
        MASKNAME = 'maskC'
        Nt,Nj,Ni = (13,270,270)

    # read ECCO tile grid file
    invariant = ncdf_invariant(os.path.join(DIRECTORY,'ECCO-GRID.nc'),
        lon=LONNAME, lat=LATNAME, depth=ZNAME, area=AREANAME, mask=MASKNAME)
    # bad value
    fill_value = -1e+10
    # model gamma and rhonil
    gamma = 9.81
    rhonil = 1029

    # compile regular expression operator for finding files for years
    regex_years = r'|'.join(rf'{y:d}' for y in range(RANGE[0],RANGE[1]+1))
    rx1 = re.compile(rf'PHIBOT([\.\_])({regex_years})(_(\d+))?.nc$')
    # find input files
    input_files = [fi for fi in os.listdir(DIRECTORY) if rx1.match(fi)]

    # output multi-annual mean
    obp_mean = {}
    # allocate for output data
    obp_mean[VARNAME] = np.ma.zeros((Nt,Nj,Ni),fill_value=fill_value)
    obp_mean[VARNAME].mask = np.logical_not(invariant['mask'][0,:,:,:]) | \
        (invariant['depth'] == 0.0)
    # create a float object for times
    obp_mean['time'] = 0.0
    # copy geolocation variables
    obp_mean['lon'] = np.copy(invariant['lon'])
    obp_mean['lat'] = np.copy(invariant['lat'])
    # counter variable for dates
    count = 0.0
    # read each input file
    for t,fi in enumerate(input_files):
        # Open netCDF4 datafile for reading
        fileID = netCDF4.Dataset(os.path.join(DIRECTORY,fi),'r')
        # copy grid variables
        for key in ('i','j','tile'):
            obp_mean[key] = fileID.variables[key][:].copy()
        # time within netCDF files is days since epoch
        TIME = fileID.variables[TIMENAME][:].copy()
        time_string = fileID.variables[TIMENAME].units
        epoch1,to_secs = gravtk.time.parse_date_string(time_string)
        # read ocean bottom pressure anomalies for each month
        for m,delta_time in enumerate(to_secs*TIME):
            # convert from ocean bottom pressure anomalies to absolute
            PHIBOT = fileID.variables[VARNAME][m,:,:,:].copy()
            obp_tile = invariant['depth']*rhonil*gamma + PHIBOT*rhonil

            # calculate Julian day by converting to MJD and adding offset
            JD = gravtk.time.convert_delta_time(delta_time,
                epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
                scale=1.0/86400.0) + 2400000.5
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
                FORMAT='tuple')
            # convert from calendar dates to year-decimal
            obp_mean['time'] += gravtk.time.convert_calendar_decimal(
                YY,MM,day=DD,hour=hh,minute=mm,second=ss)

            # global area average of each ocean bottom pressure map is removed
            # (Greatbatch correction) https://doi.org/10.1029/94JC00847
            total_area = 0.0
            total_newton = 0.0
            # for each tile
            for k in range(0, Nt):
                # Grid point areas (m^2)
                area = invariant['area'][k,:,:]
                # calculate the tile point weight in newtons
                newtons = obp_tile[k,:,:]*area
                # mask for tile
                mask = np.logical_not(obp_mean[VARNAME].mask[k,:,:])
                # finding ocean points at each lat
                if np.count_nonzero(mask):
                    indj,indi = np.nonzero(mask)
                    # total area
                    total_area += np.sum(area[indj,indi])
                    # total weight in newtons
                    total_newton += np.sum(newtons[indj,indi])
            # remove global area average of each OBP map
            ratio = (total_newton/total_area)
            obp_mean[VARNAME].data[:,:,:] += (obp_tile - ratio)
            count += 1.0

    # convert from totals to means
    indt,indj,indi = np.nonzero(~obp_mean[VARNAME].mask)
    obp_mean[VARNAME].data[indt,indj,indi] /= count
    obp_mean['time'] /= count
    # replace invalid values with fill value
    obp_mean[VARNAME].data[obp_mean[VARNAME].mask] = obp_mean[VARNAME].fill_value

    # Defining output attributes
    attributes = {}
    attributes['title'] = f'Mean_Ocean_Bottom_Pressure_from_ECCO_{MODEL}_Model'
    # dimension attributes
    attributes['i'] = {}
    attributes['i']['long_name'] = 'x-dimension of the t grid'
    attributes['i']['axis'] = 'X'
    attributes['j'] = {}
    attributes['j']['long_name'] = 'y-dimension of the t grid'
    attributes['j']['axis'] = 'Y'
    attributes['tile'] = {}
    attributes['tile']['long_name'] = 'index of llc grid tile'
    attributes[TIMENAME] = {}
    attributes[TIMENAME]['long_name'] = 'Date_in_Decimal_Years'
    attributes[TIMENAME]['units'] = 'years'
    # longitude and latitude
    attributes['lon'] = {}
    attributes['lon']['long_name'] = 'longitude'
    attributes['lon']['units'] = 'degrees_east'
    attributes['lat'] = {}
    attributes['lat']['long_name'] = 'latitude'
    attributes['lat']['units'] = 'degrees_north'
    # output ocean bottom pressure
    attributes[VARNAME] = {}
    attributes[VARNAME]['long_name'] = 'pressure_at_sea_floor'
    attributes[VARNAME]['units'] = 'Pa'

    # output to file
    FILE = f'ECCO_{MODEL}_OBP_MEAN_{RANGE[0]:4d}-{RANGE[1]:4d}.nc'
    # netcdf (.nc)
    ncdf_tile_write(obp_mean, attributes,
        FILENAME=os.path.join(DIRECTORY,FILE),
        LONNAME='lon', LATNAME='lat', TIMENAME=TIMENAME,
        VARNAME=VARNAME)
    # change the permissions mode of the output file to MODE
    os.chmod(os.path.join(DIRECTORY,FILE),MODE)

# PURPOSE: read ECCO invariant grid file
def ncdf_invariant(invariant_file,**kwargs):
    # output dictionary with invariant parameters
    invariant = {}
    # open netCDF4 file for reading
    with netCDF4.Dataset(os.path.expanduser(invariant_file),'r') as fileID:
        # extract latitude, longitude, depth, area and valid mask
        for key,val in kwargs.items():
            invariant[key] = fileID.variables[val][:].copy()
    # return the invariant parameters
    return invariant

# PURPOSE: write tiled data to a netCDF4 flie
def ncdf_tile_write(output, attributes, FILENAME=None, LONNAME=None,
    LATNAME=None, TIMENAME=None, VARNAME=None):

    # opening NetCDF file for writing
    fileID = netCDF4.Dataset(os.path.expanduser(FILENAME),'w')

    # python dictionary with NetCDF variables
    nc = {}
    # Defining the NetCDF dimensions and variables
    for key in ('i','j','tile',TIMENAME):
        fileID.createDimension(key, len(np.atleast_1d(output[key])))
        nc[key] = fileID.createVariable(key, output[key].dtype, (key,))
        # filling NetCDF variables
        nc[key][:] = np.copy(output[key])
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            setattr(nc[key],att_name,att_val)

    # Defining the NetCDF variables
    for key in (LONNAME,LATNAME,VARNAME):
        if hasattr(output[key],'fill_value'):
            nc[key] = fileID.createVariable(key, output[key].dtype,
                ('tile','j','i'), fill_value=output[key].fill_value,
                zlib=True)
        else:
            nc[key] = fileID.createVariable(key, output[key].dtype,
                ('tile','j','i'))
        # filling NetCDF variables
        nc[key][:] = np.copy(output[key])
        # Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            setattr(nc[key],att_name,att_val)
    # add attribute for date created
    fileID.date_created = datetime.datetime.now().isoformat()
    fileID.title = attributes['title']
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # Output NetCDF structure information
    logging.info(FILENAME)
    logging.info(list(fileID.variables.keys()))
    # Closing the NetCDF file
    fileID.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads monthly ECCO ocean bottom pressure
            LLC tile data and calculates multi-annual means
            """
    )
    # command line parameters
    parser.add_argument('model',
        type=str, nargs='+',
        default=['V4r4','V5alpha'], choices=['V4r4','V5alpha'],
        help='ECCO Version 4 or 5 Model')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[2003,2007],
        help='Start and end year range for mean')
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

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # for each ECCO LLC tile model
    for MODEL in args.model:
        # run program
        ecco_mean_llc_tiles(args.directory, MODEL, RANGE=args.mean,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
