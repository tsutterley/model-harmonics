#!/usr/bin/env python
u"""
era5_land_mean_monthly.py
Written by Tyler Sutterley (04/2025)
Reads monthly ERA5-land datafiles to calculate multi-annual means of
terrestrial water storage from soil moisture, snow water equivalent
and total canopy storage

snow_depth_water_equivalent (sd)
snow_cover (snowc): 0 -- 100 percent
skin_reservoir_content (src)
volumetric_soil_water_layer_1 (swvl1): 0 -- 7 cm
volumetric_soil_water_layer_2 (swvl2): 7 -- 28 cm
volumetric_soil_water_layer_3 (swvl3): 28 -- 100 cm
volumetric_soil_water_layer_4 (swvl4): 100 -- 289 cm

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -m X, --mean X: Year range for mean
    -F X, --format X: Input and output data format
        ascii
        netCDF4
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
    Written 04/2025
"""
from __future__ import print_function

import sys
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk

# ERA5-land products
variables = ['sd','snowc','src','swvl1','swvl2','swvl3','swvl4']

# PURPOSE: read variables from ERA5-Land files
def read_era5_variables(era5_land_file, **kwargs):
    # set default variables
    kwargs.setdefault('variables', variables)
    # python dictionary of output variables
    dinput = {}
    attrs = {}
    # read each variable of interest in ERA5-Land file
    era5_land_file = pathlib.Path(era5_land_file).expanduser().absolute()
    logging.debug(str(era5_land_file))
    with netCDF4.Dataset(era5_land_file, mode='r') as fileID:
        # extract geolocation variables
        dinput['latitude'] = fileID.variables['latitude'][:].copy()
        dinput['longitude'] = fileID.variables['longitude'][:].copy()
        # convert time from netCDF4 units to Julian Days
        date_string = fileID.variables['valid_time'].units
        epoch,to_secs = gravtk.time.parse_date_string(date_string)
        dinput['time'] = gravtk.time.convert_delta_time(
            to_secs*fileID.variables['valid_time'][:],epoch1=epoch,
            epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0) + 2400000.5
        # read each variable of interest in ERA5-Land file
        for key in kwargs['variables']:
            # Getting the data from each NetCDF variable of interest
            # check dimensions for expver slice
            if (fileID.variables[key].ndim == 4):
                dinput[key] = ncdf_expver(fileID, key)
            else:
                dinput[key] = np.ma.array(fileID.variables[key][:].squeeze(),
                    fill_value=fileID.variables[key]._FillValue)
            dinput[key].mask = (dinput[key].data == dinput[key].fill_value)
            # get attributes for variable
            attrs[key] = {}
            for att_name in fileID.variables[key].ncattrs():
                attrs[key][att_name] = fileID.variables[key].getncattr(att_name)
    # return the output variables
    return (dinput, attrs)

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

def era5_land_mean_monthly(base_dir,
        RANGE=None,
        DATAFORM=None,
        VERBOSE=False,
        MODE=0o775
    ):

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[VERBOSE])

    # attributes for output files
    attributes = {}
    attributes['units'] = 'm w.e.'
    attributes['longname'] = 'Equivalent_Water_Thickness'
    attributes['title'] = 'ERA5-Land Terrestrial Water Storage'
    attributes['source'] = ', '.join(variables)
    attributes['reference'] = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # directory models
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # ERA5-land products
    MODEL = 'ERA5-Land'
    ddir = base_dir.joinpath(MODEL)
    # output data file format
    suffix = dict(ascii='txt', netCDF4='nc', HDF5='H5')[DATAFORM]
    # years to run for mean
    years = np.arange(RANGE[0], RANGE[1]+1)
    # output dimensions
    nlat,nlon = (1801,3600)
    # output bad value
    fill_value = -9999.0

    # allocate for TWS and date
    tws = gravtk.spatial(fill_value=fill_value)
    tws.data = np.zeros((nlat,nlon))
    tws.mask = np.zeros((nlat,nlon),dtype=bool)
    tws.time = 0.0
    # create counter for dates
    c = 0
    # for each year within years_range
    for i,y in enumerate(years):
        # read ERA5-Land data
        model_file = f"{MODEL}-Monthly-{y:4d}.nc"
        dinput, attrs = read_era5_variables(ddir.joinpath(model_file))
        # copy the geolocation variables
        tws.lon = dinput['longitude'].copy()
        tws.lat = dinput['latitude'].copy()
        # iterate over months
        for m in range(12):
            # skin reservoir content (canopy water)
            CW = dinput['src'][m,:,:]
            # calculate snow water equivalent (SWE)
            # convert snow cover from percent to fraction
            SWE = dinput['sd'][m,:,:]*dinput['snowc'][m,:,:]/100.0
            # convert volumetric soil water content to m w.e.
            SML1 = (0.07 - 0.0)*dinput['swvl1'][m,:,:]
            SML2 = (0.28 - 0.07)*dinput['swvl2'][m,:,:]
            SML3 = (1.0 - 0.28)*dinput['swvl3'][m,:,:]
            SML4 = (2.89 - 1.0)*dinput['swvl4'][m,:,:]
            # calculate terrestrial water storage (TWS)
            tws.data[:,:] += SML1 + SML2 + SML3 + SML4 + SWE + CW
            # set the mask for invalid points
            tws.mask[:,:] |= np.isnan(CW) | (CW == attrs['src']['_FillValue'])
            tws.replace_masked()
            # convert from Julian days to calendar dates
            YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(
                dinput['time'][m], FORMAT='tuple')
            # convert from calendar dates to year-decimal
            tws.time += gravtk.time.convert_calendar_decimal(
                YY,MM,day=DD,hour=hh,minute=mm,second=ss)
            # add 1 to counter
            c += 1

    # calculate mean TWS and date
    tws_mean = tws.scale(1.0/c)
    tws_mean.time = tws.time/c

    # output to file
    model_file = f"{MODEL}-TWS-Mean-{RANGE[0]:4d}-{RANGE[1]:4d}.{suffix}"
    output_file = ddir.joinpath(model_file)

    if (DATAFORM == 'ascii'):
        # ascii (.txt)
        tws_mean.to_ascii(output_file, date=False,
            verbose=VERBOSE)
    elif (DATAFORM == 'netCDF4'):
        # netCDF4 (.nc)
        tws_mean.to_netCDF4(output_file, verbose=VERBOSE,
            date=False, **attributes)
    elif (DATAFORM == 'HDF5'):
        # HDF5 (.H5)
        tws_mean.to_HDF5(output_file, verbose=VERBOSE,
            date=False, **attributes)
    # change the permissions mode
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads ERA5-Land datafiles to calculate the
            multi-annual mean terrestrial water storage
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[2003,2007],
        help='Start and end year range for mean')
    # input and output data format (ascii, netCDF4, HDF5)
    parser.add_argument('--format','-F',
        type=str, default='netCDF4', choices=['ascii','netCDF4','HDF5'],
        help='Input and output data format')
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

    # run program
    era5_land_mean_monthly(args.directory,
        RANGE=args.mean, DATAFORM=args.format,
        VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
