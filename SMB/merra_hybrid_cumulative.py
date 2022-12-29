#!/usr/bin/env python
u"""
merra_hybrid_cumulative.py
Written by Tyler Sutterley (12/2022)
Reads MERRA-2 hybrid datafiles to calculate cumulative anomalies in
    derived surface mass balance products
MERRA-2 Hybrid model outputs provided by Brooke Medley at GSFC

CALLING SEQUENCE:
    python merra_hybrid_cumulative.py --directory <path> --region gris \
        --mean 1980 1995

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -R X, --region X: Region to calculate (gris, ais)
    -v X, --version X: Version of firn model to calculate
        v0
        v1
        v1.0
        v1.1
        v1.2
        v1.2.1
    --mean: Start and end year of mean
    -G, --gzip: netCDF4 file is locally gzip compressed
    -V, --verbose: Output information for each output file
    -M X, --mode X: Local permissions mode of the directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: add Greenland and Antarctic versions v1.2.1
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 12/2021: added GSFC MERRA-2 Hybrid Greenland v1.2
        can use variable loglevels for verbose output
    Updated 10/2021: using python logging for handling verbose output
    Updated 08/2021: output areas to file if applicable
        add verbose option to print input and output file information
        additionally output surface mass balance anomalies
    Updated 02/2021: using argparse to set parameters
        read and write for all available variables in a file
        added gzip compression option
    Written 10/2019
"""
from __future__ import print_function

import os
import re
import gzip
import time
import uuid
import logging
import netCDF4
import argparse
import numpy as np
import model_harmonics as mdlhmc

# PURPOSE: calculate cumulative anomalies in MERRA-2 hybrid
# surface mass balance variables
def merra_hybrid_cumulative(base_dir, REGION, VERSION,
    RANGE=None, GZIP=False, MODE=0o775):
    """
    Calculates cumulative anomalies of MERRA-2 hybrid
    surface mass balance products

    Parameters
    ----------
    base_dir: str
        Working data directory
    REGION: str
        MERRA-2 region to interpolate

            - ``ais``: Antarctica
            - ``gris``: Greenland
    VERSION: str
        MERRA-2 hybrid model version
    RANGE: list
        Start and end year for mean
    GZIP: bool, default False
        netCDF4 file is gzip compressed
    VERBOSE: bool, default False
        Verbose output of netCDF4 variables
    MODE: oct, default 0o775
        Permission mode of directories and files created
    """

    # MERRA-2 hybrid directory
    DIRECTORY = os.path.join(base_dir,VERSION)
    # set version parameters
    suffix = '.gz' if GZIP else ''
    if (VERSION == 'v0'):
        # input and output netCDF4 files
        args = (REGION.lower(),suffix)
        hybrid_file = 'm2_hybrid_p_minus_e_melt_{0}.nc{1}'.format(*args)
        output_file = 'm2_hybrid_cumul_{0}.nc{1}'.format(*args)
        # names of variables to read
        VARIABLES = ('p_minus_e','melt')
        AREA = None
        anomaly_flag = '_anomaly'
    elif VERSION in ('v1','v1.0'):
        # input and output netCDF4 files
        MAJOR_VERSION = re.match(r'((v\d+)(\.\d+)?)$',VERSION).group(2)
        args = (MAJOR_VERSION,REGION.lower(),suffix)
        hybrid_file = 'gsfc_fdm_smb_{0}_{1}.nc{2}'.format(*args)
        output_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
        # names of variables to read
        VARIABLES = ('runoff','rainfall','snowfall_minus_sublimation','SMB')
        AREA = None
        # flag to append to output netCDF4 variables
        anomaly_flag = '_anomaly'
    else:
        # input and output netCDF4 files
        FILE_VERSION = VERSION.replace('.','_')
        args = (FILE_VERSION,REGION.lower(),suffix)
        firn_height_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
        hybrid_file = 'gsfc_fdm_smb_{0}_{1}.nc{2}'.format(*args)
        output_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
        # names of variables to read
        VARIABLES = ('Me','Ra','Ru','Sn-Ev','SMB')
        AREA = 'iArea'
        # flag to append to output netCDF4 variables
        anomaly_flag = '_a'

    # Open the MERRA-2 Hybrid NetCDF file for reading
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(os.path.join(DIRECTORY,hybrid_file),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(os.path.join(DIRECTORY,hybrid_file), 'r')

    # Output NetCDF file information
    logging.info(os.path.join(DIRECTORY,hybrid_file))
    logging.info(list(fileID.variables.keys()))

    # Get data and attribute from each netCDF variable
    fd = {}
    attrs = {}
    # input time (year-decimal)
    fd['time'] = fileID.variables['time'][:].copy()
    # extract areas from SMB or firn height file
    if AREA and (AREA in fileID.variables):
        fd[AREA] = fileID.variables[AREA][:].copy()
        # get each attribute for area variable if applicable
        attrs[AREA] = {}
        for att_name in ['units','long_name','standard_name']:
            if hasattr(fileID.variables[AREA],att_name):
                attrs[AREA][att_name]=fileID.variables[AREA].getncattr(att_name)
    elif AREA:
        # Open the MERRA-2 Hybrid firn height file for reading
        if GZIP:
            # read as in-memory (diskless) netCDF4 dataset
            with gzip.open(os.path.join(DIRECTORY,firn_height_file),'r') as f:
                fid1 = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
        else:
            # read netCDF4 dataset
            fid1 = netCDF4.Dataset(os.path.join(DIRECTORY,firn_height_file), 'r')
        # copy area from firn height file
        fd[AREA] = fid1.variables[AREA][:].copy()
        # get each attribute for area variable if applicable
        attrs[AREA] = {}
        for att_name in ['units','long_name','standard_name']:
            if hasattr(fid1.variables[AREA],att_name):
                attrs[AREA][att_name]=fid1.variables[AREA].getncattr(att_name)
        # close the firn height file
        fid1.close()

    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    if (np.ndim(fileID.variables['x'][:]) == 2):
        xg = fileID.variables['x'][:].copy()
        yg = fileID.variables['y'][:].copy()
        fd['x'],fd['y'] = (xg[:,0],yg[0,:])
    else:
        fd['x'] = fileID.variables['x'][:].copy()
        fd['y'] = fileID.variables['y'][:].copy()
        xg,yg = np.meshgrid(fd['x'],fd['y'],indexing='ij')
    # time is year decimal at time step 5 days
    time_step = 5.0/365.25
    # calculate mean period for MERRA-2
    tt, = np.nonzero((fd['time'] >= RANGE[0]) & (fd['time'] < (RANGE[1]+1)))

    # for each variable
    for v in VARIABLES:
        # copy data and remove singleton dimensions
        DATA = np.ma.array(fileID.variables[v][:]).squeeze()
        # invalid data value
        DATA.fill_value = np.float64(fileID.variables[v]._FillValue)
        # set masks
        DATA.mask = (DATA.data == DATA.fill_value)
        # get each attribute for variable if applicable
        attrs[v] = {}
        for att_name in ['units','long_name','standard_name','comment']:
            if hasattr(fileID.variables[v],att_name):
                attrs[v][att_name] = fileID.variables[v].getncattr(att_name)
        # input shape of MERRA-2 Hybrid firn data
        nt,nx,ny = np.shape(DATA)

        # cumulative mass anomalies calculated by removing mean balance flux
        # mean of data for variable (converted from yearly rate)
        MEAN = np.mean(DATA.data[tt,:,:]*time_step, axis=0)
        # indices of specified ice mask at the first slice
        i,j = np.nonzero(~DATA.mask[0,:,:])
        valid_count = np.count_nonzero(~DATA.mask[0,:,:])
        # allocate for output variable
        fd[v] = np.ma.zeros((nt,nx,ny),fill_value=DATA.fill_value)
        fd[v].mask = (DATA.mask | np.isnan(DATA.data))
        CUMULATIVE = np.zeros((valid_count))
        # calculate output cumulative anomalies for variable
        for t in range(nt):
            # convert mass flux from yearly rate and
            # calculate cumulative anomalies at time t
            CUMULATIVE += (DATA.data[t,i,j]*time_step - MEAN[i,j])
            fd[v].data[t,i,j] = CUMULATIVE.copy()
        # replace masked values with fill value
        fd[v].data[fd[v].mask] = fd[v].fill_value
    # close the NetCDF files
    fileID.close()

    # Output NetCDF filename
    logging.info(os.path.join(DIRECTORY,output_file))

    # output MERRA-2 data file with cumulative data
    if GZIP:
        # open virtual file object for output
        fileID = netCDF4.Dataset(uuid.uuid4().hex,'w',memory=True,
            format='NETCDF4')
    else:
        # opening NetCDF file for writing
        fileID = netCDF4.Dataset(os.path.join(DIRECTORY,output_file),'w',
            format="NETCDF4")

    # Defining the NetCDF dimensions
    fileID.createDimension('x', nx)
    fileID.createDimension('y', ny)
    fileID.createDimension('time', nt)

    # python dictionary with netCDF4 variables
    nc = {}
    # defining the NetCDF variables
    nc['x'] = fileID.createVariable('x', fd['x'].dtype, ('x',))
    nc['y'] = fileID.createVariable('y', fd['y'].dtype, ('y',))
    nc['time'] = fileID.createVariable('time', fd['time'].dtype, ('time',))
    # output area variable
    if AREA:
        nc[AREA] = fileID.createVariable(AREA, fd[AREA].dtype, ('x','y',),
            fill_value=fd[AREA].fill_value, zlib=True)
    # for each output variable
    for v in VARIABLES:
        # append anomaly flag
        var = f'{v}{anomaly_flag}'
        nc[v] = fileID.createVariable(var, fd[v].dtype, ('time','x','y',),
            fill_value=fd[v].fill_value, zlib=True)

    # filling NetCDF variables
    for key,val in fd.items():
        nc[key][:] = val.copy()

    # create variable and attributes for projection
    if REGION in ('gris',):
        crs = fileID.createVariable('Polar_Stereographic',np.byte,())
        crs.standard_name = 'Polar_Stereographic'
        crs.grid_mapping_name = 'polar_stereographic'
        crs.straight_vertical_longitude_from_pole = -45.0
        crs.latitude_of_projection_origin = 90.0
        crs.standard_parallel = 70.0
        crs.scale_factor_at_projection_origin = 1.
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        crs.semi_major_axis = 6378.137
        crs.semi_minor_axis = 6356.752
        crs.inverse_flattening = 298.257223563
        crs.spatial_epsg = '3413'
    elif REGION in ('ais',):
        crs = fileID.createVariable('Polar_Stereographic',np.byte,())
        crs.standard_name = 'Polar_Stereographic'
        crs.grid_mapping_name = 'polar_stereographic'
        crs.straight_vertical_longitude_from_pole = 0.0
        crs.latitude_of_projection_origin = -90.0
        crs.standard_parallel = -71.0
        crs.scale_factor_at_projection_origin = 1.
        crs.false_easting = 0.0
        crs.false_northing = 0.0
        crs.semi_major_axis = 6378.137
        crs.semi_minor_axis = 6356.752
        crs.inverse_flattening = 298.257223563
        crs.spatial_epsg = '3031'

    # Defining attributes for x and y coordinates
    nc['x'].long_name = 'Easting'
    nc['x'].standard_name = 'projection_x_coordinate'
    nc['x'].grid_mapping = 'Polar_Stereographic'
    nc['x'].units = 'meters'
    nc['y'].long_name = 'Northing'
    nc['y'].standard_name = 'projection_y_coordinate'
    nc['y'].grid_mapping = 'Polar_Stereographic'
    nc['y'].units = 'meters'
    # defining attributes for area variable
    if AREA:
        # set area variable attributes
        for att_name,att_val in attrs[AREA].items():
            nc[AREA].setncattr(att_name,att_val)
        # set grid mapping attribute
        nc[AREA].setncattr('grid_mapping','Polar_Stereographic')
    # Defining attributes for variables
    for v in VARIABLES:
        # set variable attributes
        for att_name,att_val in attrs[v].items():
            nc[v].setncattr(att_name,att_val.replace(' per year',''))
        # set grid mapping attribute
        nc[v].setncattr('grid_mapping','Polar_Stereographic')
    # Defining attributes for date
    nc['time'].long_name = 'time, 5-daily resolution'
    nc['time'].units = 'decimal years, 5-daily resolution'
    # global attributes of NetCDF file
    fileID.title = (f'Cumulative anomalies in GSFC-FDM{VERSION} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}')
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
    fileID.source = f'version {VERSION}'
    fileID.references = ("Medley, B., Neumann, T. A., Zwally, H. J., "
        "Smith, B. E., and Stevens, C. M.: Simulations of Firn Processes "
        "over the Greenland and Antarctic Ice Sheets: 1980--2021, "
        "The Cryosphere, https://doi.org/10.5194/tc-2020-266, 2022.")
    fileID.institution = "NASA Goddard Space Flight Center (GSFC)"
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    # Output NetCDF file information
    logging.info(list(fileID.variables.keys()))
    # Closing the NetCDF file and getting the buffer object
    nc_buffer = fileID.close()

    # write MERRA-2 data file to gzipped file
    if GZIP:
        # copy bytes to file
        with gzip.open(os.path.join(DIRECTORY,output_file), 'wb') as f:
            f.write(nc_buffer)

    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,output_file), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads MERRA-2 Hybrid datafiles to
            calculate cumulative anomalies in surface
            mass balance products
            """
    )
    # command line parameters
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # region of firn model
    parser.add_argument('--region','-R',
        type=str, default='gris', choices=['gris','ais'],
        help='Region of firn model to calculate')
    # version of firn model
    versions = ['v0','v1','v1.0','v1.1','v1.2','v1.2.1']
    parser.add_argument('--version','-v',
        type=str, default='v1.2.1', choices=versions,
        help='Version of firn model to calculate')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    # netCDF4 files are gzip compressed
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='netCDF4 file is locally gzip compressed')
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

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    merra_hybrid_cumulative(args.directory, args.region, args.version,
        RANGE=args.mean, GZIP=args.gzip, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
