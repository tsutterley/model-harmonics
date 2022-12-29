#!/usr/bin/env python
u"""
racmo_smb_cumulative.py
Written by Tyler Sutterley (12/2022)
Reads RACMO datafiles to calculate cumulative anomalies in derived surface
    mass balance products

CALLING SEQUENCE:
    python racmo_smb_cumulative.py --product smb --verbose <path_to_racmo_file>

COMMAND LINE OPTIONS:
    -P X, --product X: RACMO SMB product to calculate
    --mean: Start and end year of mean
    -G, --gzip: netCDF4 file is locally gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Local permissions mode of the directories and files

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
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 04/2022: deprecation fixes for regular expressions
    Updated 12/2021: can use variable loglevels for verbose output
    Updated 11/2021: complete rewrite of program
        dropped old RACMO ascii file read portions
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: concatenated numpy arrays as input to nearest neighbors
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018: can use getopt to set the working directory
    Updated 01/2018: using scikit-learn nearest neighbors for quickly finding
        mapping between input model grid and output global grid
    Updated 10/2017: new way of selecting between model versions
    Updated 09/2017: added RACMO2.3 XANT27 1979-2016
    Updated 02/2017: base_dir as input to racmo functions
    Updated 07/2016: added RACMO2.3 1958-2015 model, using netCDF4-python
    Updated 06/2016: using __future__ print function
    Updated 12/2015: RACMO2.3 XANT27 run for 1979-2014
    Updated 06/2015: added new PRODUCTS including derived Rainfall product
        saving indices for converting from RACMO grid to global grid
    Updated 04/2015: added HDF5 output option
    Updated 03/2015: updated for new test_ant27 program
    Updated 06/2014: can run from command line and parameter file
    Updated 05/2014: code generalizations and updates
    Updated 02/2014: more general code updates
    Updated 01/2014: create date files. add new RACMO 2.3 model
    Written 10/2011
"""
from __future__ import print_function

import sys
import os
import re
import gzip
import uuid
import time
import logging
import netCDF4
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: read and cumulative RACMO SMB SMB estimates
def racmo_smb_cumulative(model_file, VARIABLE,
    RANGE=None,
    GZIP=False,
    MODE=0o775):

    # RACMO SMB directory
    DIRECTORY = os.path.dirname(model_file)
    # try to extract region and version from filename
    R1 = re.compile(r'[XF]?(ANT27|GRN11|GRN055|PEN055|ASE055)',re.VERBOSE)
    R2 = re.compile(r'(RACMO\d+(\.\d+)?(p\d+)?)',re.VERBOSE)
    REGION = R1.search(os.path.basename(model_file)).group(0)
    VERSION = R2.search(os.path.basename(model_file)).group(0)
    # RACMO products
    racmo_products = {}
    racmo_products['precip'] = 'Precipitation'
    racmo_products['rainfall'] = 'Rainfall'
    racmo_products['refreeze'] = 'Meltwater Refreeze'
    racmo_products['runoff'] = 'Meltwater Runoff'
    racmo_products['smb'] = 'Surface Mass Balance'
    racmo_products['sndiv'] = 'Snow Drift Erosion'
    racmo_products['snowfall'] = 'Snowfall'
    racmo_products['snowmelt'] = 'Snowmelt'
    racmo_products['subl'] = 'Sublimation'

    # Open the RACMO SMB NetCDF file for reading
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.expanduser(model_file),'r') as f:
            f_in = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # Open the RACMO NetCDF file for reading
        f_in = netCDF4.Dataset(os.path.expanduser(model_file), 'r')

    # Output NetCDF file information
    logging.info(os.path.expanduser(model_file))
    logging.info(list(f_in.variables.keys()))

    # Get data and attributes for each netCDF variable
    fd,attrs = ({},{})
    attributes_list = ['axis','calendar','description','grid_mapping',
        'long_name','standard_name','units','_FillValue']
    for key in [VARIABLE,'lon','lat','rlon','rlat','time']:
        # remove singleton dimensions
        fd[key] = np.squeeze(f_in.variables[key][:].copy())
        # get applicable attributes for variable
        attrs[key] = {}
        for att_name in attributes_list:
            # try getting the attribute
            try:
                ncattr, = [s for s in f_in[key].ncattrs() if
                    re.match(att_name,s,re.I)]
                attrs[key][att_name] = f_in[key].getncattr(ncattr)
            except (ValueError,AttributeError):
                pass
            else:
                if isinstance(attrs[key][att_name],str):
                    attrs[key][att_name] = attrs[key][att_name].strip()

    # parse date string within netCDF4 file
    date_string = attrs['time']['units']
    epoch1,to_secs = gravtk.time.parse_date_string(date_string)
    # calculate Julian day by converting to MJD and adding offset
    JD = gravtk.time.convert_delta_time(fd['time']*to_secs,
        epoch1=epoch1, epoch2=(1858,11,17,0,0,0),
        scale=1.0/86400.0) + 2400000.5
    # convert from Julian days to calendar dates
    YY,MM,DD,hh,mm,ss = gravtk.time.convert_julian(JD,
        FORMAT='tuple')
    # convert from calendar dates to year-decimal
    TIME = gravtk.time.convert_calendar_decimal(YY,MM,
        day=DD,hour=hh,minute=mm,second=ss)

    # copy data to masked array
    DATA = np.ma.array(fd[VARIABLE].copy())
    # invalid data value
    DATA.fill_value = np.float64(attrs[VARIABLE]['_FillValue'])
    # set masks
    DATA.mask = (DATA.data == DATA.fill_value)
    # input shape of RACMO data
    nt,ny,nx = np.shape(DATA)

    # calculate mean period for RACMO
    tt, = np.nonzero((TIME >= RANGE[0]) & (TIME < (RANGE[1]+1)))
    # cumulative mass anomalies calculated by removing mean balance flux
    MEAN = np.mean(DATA.data[tt,:,:], axis=0)
    # indices of specified ice mask at the first slice
    i,j = np.nonzero(~DATA.mask[0,:,:])
    valid_count = np.count_nonzero(~DATA.mask[0,:,:])
    # allocate for output variable
    fd[VARIABLE] = np.ma.zeros((nt,ny,nx),fill_value=DATA.fill_value)
    fd[VARIABLE].mask = (DATA.mask | np.isnan(DATA.data))
    CUMULATIVE = np.zeros((valid_count))
    # calculate output cumulative anomalies for variable
    for t in range(nt):
        # convert mass flux from yearly rate and
        # calculate cumulative anomalies at time t
        CUMULATIVE += (DATA.data[t,i,j] - MEAN[i,j])
        fd[VARIABLE].data[t,i,j] = CUMULATIVE.copy()
    # replace masked values with fill value
    fd[VARIABLE].data[fd[VARIABLE].mask] = fd[VARIABLE].fill_value

    # Output NetCDF filename
    FILE = f'{VERSION}_{REGION}_{VARIABLE.upper()}_cumul.nc'
    logging.info(os.path.join(DIRECTORY,FILE))

    # output MERRA-2 data file with cumulative data
    if GZIP:
        # open virtual file object for output
        f_out = netCDF4.Dataset(uuid.uuid4().hex,'w',memory=True,
            format='NETCDF4')
    else:
        # opening NetCDF file for writing
        f_out = netCDF4.Dataset(os.path.join(DIRECTORY,FILE),'w',
            format="NETCDF4")

    # python dictionary with netCDF4 variables
    nc = {}
    # defining the NetCDF dimensions
    for key in ['rlon','rlat','time']:
        f_out.createDimension(key, len(fd[key]))
        nc[key] = f_out.createVariable(key, fd[key].dtype, (key,))
    # for each geolocation variable
    for key in ['lon','lat']:
        nc[key] = f_out.createVariable(key, fd[key].dtype,
            ('rlat','rlon',), zlib=True)
    # output variable
    nc[VARIABLE] = f_out.createVariable(VARIABLE, fd[VARIABLE].dtype,
        ('time','rlat','rlon',), fill_value=DATA.fill_value, zlib=True)

    # copy variable and attributes for projection
    crs = f_out.createVariable('rotated_pole',np.byte,())
    for att_name in f_in['rotated_pole'].ncattrs():
        att_val = f_in['rotated_pole'].getncattr(att_name)
        crs.setncattr(att_name,att_val)

    # filling NetCDF variables
    for key,val in fd.items():
        nc[key][:] = val.copy()
        # for each variable attribute
        for att_name,att_val in attrs[key].items():
            if att_name not in ('_FillValue',):
                nc[key].setncattr(att_name,att_val)

    # global attributes of NetCDF file
    for att_name in ['comment','Domain','Experiment','source','title']:
        try:
            ncattr, = [s for s in f_in.variables[key].ncattrs()
                if re.match(att_name,s,re.I)]
            attribute = f_in.getncattr(ncattr)
        except (ValueError,AttributeError):
            pass
        else:
            f_out.setncattr(att_name,attribute)
    # output attribute for mean
    f_out.description = (f'Cumulative anomalies in {VERSION} {REGION} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}')
    # add software information
    f_out.software_reference = mdlhmc.version.project_name
    f_out.software_version = mdlhmc.version.full_version
    f_out.reference = f'Output from {os.path.basename(sys.argv[0])}'
    # date created
    f_out.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output NetCDF file information
    logging.info(list(f_out.variables.keys()))

    # Closing the NetCDF file and getting the buffer object
    f_in.close()
    nc_buffer = f_out.close()

    # write RACMO data file to gzipped file
    if GZIP:
        # copy bytes to file
        with gzip.open(os.path.join(DIRECTORY,FILE), 'wb') as f:
            f.write(nc_buffer)

    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,FILE), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates cumulative anomalies of RACMO
            surface mass balance products
            """
    )
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='RACMO SMB file to run')
    # products from SMB model
    choices = ['precip','rainfall','refreeze','runoff','smb',
        'sndiv','snowfall','snowmelt','subl']
    parser.add_argument('--product','-P',
        type=str, metavar='PRODUCT', default='smb', choices=choices,
        help='RACMO SMB product to calculate')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    # netCDF4 files are gzip compressed
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='netCDF4 file is locally gzip compressed')
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

    # run program
    racmo_smb_cumulative(args.infile, args.product,
        RANGE=args.mean,
        GZIP=args.gzip,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()