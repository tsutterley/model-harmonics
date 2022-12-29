#!/usr/bin/env python
u"""
gemb_smb_cumulative.py
Written by Tyler Sutterley (12/2022)
Calculates cumulative anomalies of GEMB surface mass balance products

CALLING SEQUENCE:
    python gemb_smb_cumulative.py --mean 1980 1995 <path_to_gemb_file>

COMMAND LINE OPTIONS:
    --mean: Start and end year of mean
    -f X, --fill-value X: set fill_value for input spatial fields
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
    Written 10/2022
"""
from __future__ import print_function

import os
import re
import time
import logging
import netCDF4
import argparse
import numpy as np
import model_harmonics as mdlhmc

# PURPOSE: calculate cumulative anomalies in GEMB
# surface mass balance variables
def gemb_smb_cumulative(model_file,
    RANGE=None,
    FILL_VALUE=np.nan,
    MODE=0o775):
    """
    Calculates cumulative anomalies of GEMB
    surface mass balance products

    Parameters
    ----------
    model_file: str
        input GEMB file
    RANGE: list
        Start and end year for mean
    FILL_VALUE: float, default np.nan
        Output invalid value
    MODE: oct, default 0o775
        Permission mode of directories and files created
    """

    # GEMB directory
    DIRECTORY = os.path.dirname(model_file)
    # regular expression pattern for extracting parameters
    pattern = r'GEMB_(Greenland|Antarctica)_SMB_\d{4}_\d{4}_mesh_\d+km_(v.*?).nc$'
    region, version = re.findall(pattern, model_file).pop()
    output_file = f'GEMB_{region}_SMB_cumul_{version}.nc'

    # Open the GEMB NetCDF file for reading
    fileID = netCDF4.Dataset(os.path.expanduser(model_file), 'r')

    # Output NetCDF file information
    logging.info(os.path.expanduser(model_file))
    logging.info(list(fileID.variables.keys()))

    # Get data and attribute from each netCDF variable
    fd = {}
    attrs = {}
    # input time (year-decimal)
    fd['time'] = fileID.variables['time'][:].copy()
    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    fd['x'] = fileID.variables['x'][:].copy()
    fd['y'] = fileID.variables['y'][:].copy()
    # calculate mean period for GEMB
    tt, = np.nonzero((fd['time'] >= RANGE[0]) & (fd['time'] < (RANGE[1]+1)))

    # copy data and remove singleton dimensions
    centered_SMB = fileID.variables['centered_SMB'][:].copy()
    accum_SMB = fileID.variables['accum_SMB'][:].copy()
    # get each attribute for variable if applicable
    for v in ['accum_SMB','centered_SMB']:
        attrs[v] = {}
        for att_name in ['units','long_name','standard_name','comment']:
            if hasattr(fileID.variables[v],att_name):
                attrs[v][att_name] = fileID.variables[v].getncattr(att_name)
    # edit cumulative SMB attributes
    attrs['accum_SMB']['standard_name'] = 'accumulated surface mass balance height'
    # input shape of GEMB SMB data
    nt,ny,nx = np.shape(accum_SMB)
    # get root attributes
    institution = fileID.getncattr('institution')
    revision = fileID.getncattr('revision')
    # close the NetCDF files
    fileID.close()

    # indices of specified ice mask
    i,j = np.nonzero(np.isfinite(centered_SMB))
    valid_count = np.count_nonzero(np.isfinite(centered_SMB))
    # calculate original monthly SMB data from anomalies
    SMB = np.ma.zeros((nt,ny,nx), fill_value=FILL_VALUE)
    SMB.mask = np.ones((nt,ny,nx), dtype=bool)
    SMB_anomaly = np.zeros((valid_count))
    # for each date
    for t in range(nt):
        SMB.data[t,i,j] = accum_SMB.data[t,i,j] - SMB_anomaly
        # set masks
        SMB.mask[t,i,j] = np.isnan(accum_SMB.data[t,i,j])
        # update SMB anomaly variable
        SMB_anomaly[:] = np.copy(accum_SMB[t,i,j])
    # convert invalid values to fill value
    SMB.data[SMB.mask] = SMB.fill_value

    # cumulative mass anomalies calculated by removing mean balance flux
    MEAN = np.mean(SMB[tt,:,:], axis=0)
    # allocate for output variable
    fd['accum_SMB'] = np.ma.zeros((nt,ny,nx), fill_value=FILL_VALUE)
    fd['accum_SMB'].mask = (SMB.mask | np.isnan(SMB.data))
    CUMULATIVE = np.zeros((valid_count))
    # calculate output cumulative anomalies for variable
    for t in range(nt):
        # calculate cumulative anomalies at time t
        CUMULATIVE += (SMB.data[t,i,j] - MEAN[i,j])
        fd['accum_SMB'].data[t,i,j] = CUMULATIVE.copy()
    # replace masked values with fill value
    fd['accum_SMB'].data[fd['accum_SMB'].mask] = fd['accum_SMB'].fill_value

    # Output NetCDF filename
    logging.info(os.path.join(DIRECTORY,output_file))

    # output GEMB data file with cumulative data
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
    nc['accum_SMB'] = fileID.createVariable('accum_SMB', fd['accum_SMB'].dtype,
        ('time','y','x',), fill_value=fd['accum_SMB'].fill_value, zlib=True)

    # filling NetCDF variables
    for key,val in fd.items():
        nc[key][:] = val.copy()

    # create variable and attributes for projection
    if region in ('Greenland',):
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
    elif region in ('Antarctica',):
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
    # set variable attributes
    for att_name,att_val in attrs['accum_SMB'].items():
        nc['accum_SMB'].setncattr(att_name,att_val)
    # set grid mapping attribute
    nc['accum_SMB'].setncattr('grid_mapping','Polar_Stereographic')
    # Defining attributes for date
    nc['time'].long_name = 'time'
    nc['time'].units = 'decimal years'
    # global attributes of NetCDF file
    fileID.title = (f'Cumulative anomalies in GEMB{version} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}')
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())
    fileID.source = f'version {version}'
    fileID.authors = "Nicole-Jeanne Schlegel & Alex Gardner"
    fileID.reference = "https://doi.org/10.5281/zenodo.7199528"
    fileID.institution = institution
    fileID.revision = revision
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    # Output NetCDF file information
    logging.info(list(fileID.variables.keys()))
    # Closing the NetCDF file and getting the buffer object
    fileID.close()

    # change the permissions mode
    os.chmod(os.path.join(DIRECTORY,output_file), MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads GEMB datafiles to calculate monthly
            cumulative anomalies in surface mass balance products
            """
    )
    # command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='GEMB file to run')
    # start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[1980,1995],
        help='Start and end year range for mean')
    # fill value for ascii
    parser.add_argument('--fill-value','-f',
        type=float, default=np.nan,
        help='Output invalid value')
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    # return the parser
    return parser

# Main program that calls gemb_smb_cumulative()
def main():
    parser = arguments()
    args = parser.parse_args()

    # create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    # run program
    gemb_smb_cumulative(args.infile,
        RANGE=args.mean,
        FILL_VALUE=args.fill_value,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
