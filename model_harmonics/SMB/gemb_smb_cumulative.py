#!/usr/bin/env python
"""
gemb_smb_cumulative.py
Written by Tyler Sutterley (10/2024)
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
    Updated 10/2024: output centered SMB in addition to cumulative
        updated regular expression for new filenames
    Updated 03/2023: regular expression pattern can find if periphery
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Written 10/2022
"""

from __future__ import print_function

import re
import time
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import model_harmonics as mdlhmc
from gravity_toolkit.utilities import build_logger

# dictionary defining the spatial reference system for each region
crs = dict(Greenland={}, Antarctica={})
# Greenland (EPSG:3413)
crs['Greenland']['standard_name'] = 'Polar_Stereographic'
crs['Greenland']['grid_mapping_name'] = 'polar_stereographic'
crs['Greenland']['straight_vertical_longitude_from_pole'] = -45.0
crs['Greenland']['latitude_of_projection_origin'] = 90.0
crs['Greenland']['standard_parallel'] = 70.0
crs['Greenland']['scale_factor_at_projection_origin'] = 1.0
crs['Greenland']['false_easting'] = 0.0
crs['Greenland']['false_northing'] = 0.0
crs['Greenland']['semi_major_axis'] = 6378.137
crs['Greenland']['semi_minor_axis'] = 6356.752
crs['Greenland']['inverse_flattening'] = 298.257223563
crs['Greenland']['spatial_epsg'] = '3413'
# Antarctica (EPSG:3031)
crs['Antarctica']['standard_name'] = 'Polar_Stereographic'
crs['Antarctica']['grid_mapping_name'] = 'polar_stereographic'
crs['Antarctica']['straight_vertical_longitude_from_pole'] = 0.0
crs['Antarctica']['latitude_of_projection_origin'] = -90.0
crs['Antarctica']['standard_parallel'] = -71.0
crs['Antarctica']['scale_factor_at_projection_origin'] = 1.0
crs['Antarctica']['false_easting'] = 0.0
crs['Antarctica']['false_northing'] = 0.0
crs['Antarctica']['semi_major_axis'] = 6378.137
crs['Antarctica']['semi_minor_axis'] = 6356.752
crs['Antarctica']['inverse_flattening'] = 298.257223563
crs['Antarctica']['spatial_epsg'] = '3031'


# PURPOSE: calculate cumulative anomalies in GEMB
# surface mass balance variables
def gemb_smb_cumulative(model_file, RANGE=None, FILL_VALUE=np.nan, MODE=0o775):
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
    # get logger
    logger = logging.getLogger(__name__)

    # regular expression pattern for extracting parameters
    pattern = (
        r'GEMB_(Greenland|Antarctica)(_and_Periphery)?_'
        r'SMB_\d{4}_\d{4}_(.*?)mesh_\d+km_(v.*?).nc$'
    )
    model_file = pathlib.Path(model_file).expanduser().absolute()
    region, aux1, aux2, version = re.findall(pattern, model_file.name).pop()
    # output cumulative file
    cumulative_file = f'GEMB_{region}{aux1}_SMB_cumul_{version}.nc'
    output_file = model_file.with_name(cumulative_file)

    # Open the GEMB NetCDF file for reading
    fileID = netCDF4.Dataset(model_file, mode='r')
    # get root attributes
    institution = fileID.getncattr('institution')
    revision = fileID.getncattr('revision')

    # Output NetCDF file information
    logger.info(str(model_file))
    logger.info(list(fileID.variables.keys()))

    # dictionary describing the output netCDF4 structure
    struct = dict(
        dimensions=('time', 'y', 'x'),
        variables={
            'accum_SMB': ('time', 'y', 'x'),
            'centered_SMB': ('y', 'x'),
        },
    )
    # dictionaries with data and attributes for output NetCDF file
    fd = {}
    attrs = dict(ROOT={})
    attributes_list = ['units', 'long_name', 'standard_name', 'comment']
    # global attributes of NetCDF file
    attrs['ROOT']['title'] = (
        f'Cumulative anomalies in GEMB{version} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}'
    )
    attrs['ROOT']['source'] = f'version {version}'
    attrs['ROOT']['authors'] = 'Nicole-Jeanne Schlegel & Alex Gardner'
    attrs['ROOT']['reference'] = 'https://doi.org/10.5281/zenodo.7199528'
    attrs['ROOT']['institution'] = institution
    attrs['ROOT']['revision'] = revision
    # list of input files for provenance
    lineage = []
    lineage.append(model_file.name)

    # create variable and attributes for projection
    fd['Polar_Stereographic'] = np.byte()
    # add projection attributes to dictionary
    attrs['Polar_Stereographic'] = crs[region]

    # input time (year-decimal)
    fd['time'] = fileID.variables['time'][:].copy()
    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    fd['x'] = fileID.variables['x'][:].copy()
    fd['y'] = fileID.variables['y'][:].copy()
    # Defining attributes for x and y coordinates
    attrs['x'] = dict(
        long_name='Easting',
        standard_name='projection_x_coordinate',
        grid_mapping='Polar_Stereographic',
        units='meters',
    )
    attrs['y'] = dict(
        long_name='Northing',
        standard_name='projection_y_coordinate',
        grid_mapping='Polar_Stereographic',
        units='meters',
    )

    # copy data and remove singleton dimensions
    centered_SMB = fileID.variables['centered_SMB'][:].copy()
    accum_SMB = fileID.variables['accum_SMB'][:].copy()
    # get each attribute for variable if applicable
    for v in ['accum_SMB', 'centered_SMB']:
        # set variable attributes
        attrs[v] = ncdf_attributes(fileID.variables[v], attributes_list)
        # set grid mapping attribute
        attrs[v]['grid_mapping'] = 'Polar_Stereographic'
    # edit cumulative SMB attributes
    attrs['centered_SMB']['standard_name'] = (
        'average surface mass balance height'
    )
    attrs['accum_SMB']['standard_name'] = (
        'accumulated surface mass balance height'
    )
    # input shape of GEMB SMB data
    nt, ny, nx = np.shape(accum_SMB)
    # close the NetCDF files
    fileID.close()

    # indices of specified ice mask
    i, j = np.nonzero(np.isfinite(centered_SMB))
    valid_count = np.count_nonzero(np.isfinite(centered_SMB))
    # calculate original monthly SMB data from anomalies
    SMB = np.ma.zeros((nt, ny, nx), fill_value=FILL_VALUE)
    SMB.mask = np.ones((nt, ny, nx), dtype=bool)
    SMB_anomaly = np.zeros((valid_count))
    # for each date
    for t in range(nt):
        SMB.data[t, i, j] = accum_SMB.data[t, i, j] - SMB_anomaly
        # set masks
        SMB.mask[t, i, j] = np.isnan(accum_SMB.data[t, i, j])
        # update SMB anomaly variable
        SMB_anomaly[:] = np.copy(accum_SMB[t, i, j])
    # convert invalid values to fill value
    SMB.data[SMB.mask] = SMB.fill_value

    # calculate mean period for GEMB
    (tt,) = np.nonzero((fd['time'] >= RANGE[0]) & (fd['time'] < (RANGE[1] + 1)))
    # cumulative mass anomalies calculated by removing mean balance flux
    MEAN = np.mean(SMB[tt, :, :], axis=0)
    fd['centered_SMB'] = np.ma.array(MEAN, fill_value=FILL_VALUE)
    fd['centered_SMB'].mask = np.isnan(MEAN)
    # allocate for output variable
    fd['accum_SMB'] = np.ma.zeros((nt, ny, nx), fill_value=FILL_VALUE)
    fd['accum_SMB'].mask = SMB.mask | np.isnan(SMB.data)
    CUMULATIVE = np.zeros((valid_count))
    # calculate output cumulative anomalies for variable
    for t in range(nt):
        # calculate cumulative anomalies at time t
        CUMULATIVE += SMB.data[t, i, j] - MEAN[i, j]
        fd['accum_SMB'].data[t, i, j] = CUMULATIVE.copy()
    # replace masked values with fill value
    fd['accum_SMB'].data[fd['accum_SMB'].mask] = fd['accum_SMB'].fill_value

    # Output NetCDF filename
    logger.info(str(output_file))
    # output GEMB data file with cumulative data
    attrs['ROOT']['lineage'] = lineage
    mdlhmc.spatial.to_netCDF4(output_file, fd, attrs, struct=struct)
    # change the permissions mode
    output_file.chmod(mode=MODE)


# PURPOSE: get attributes for netCDF4 files and variable
def ncdf_attributes(nc, attributes_list):
    # get logger
    logger = logging.getLogger(__name__)
    # output dictionary of attributes
    attributes = {}
    # for each attribute to try to get
    for att_name in attributes_list:
        rx = re.compile(att_name, re.IGNORECASE)
        try:
            # use case-insensitive regex to find attribute name
            (ncattr,) = [s for s in nc.ncattrs() if rx.match(s)]
            att_val = nc.getncattr(ncattr)
        except Exception as exc:
            ncvar = getattr(nc, 'name', 'ROOT')
            logger.debug(f'Attribute {att_name} not found in {ncvar}')
            continue
        # strip whitespace for string attributes
        if isinstance(att_val, str):
            att_val = att_val.strip()
        # add attribute to dictionary
        attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Reads GEMB datafiles to calculate monthly
            cumulative anomalies in surface mass balance products
            """
    )
    # command line parameters
    parser.add_argument('infile', type=pathlib.Path, help='GEMB file to run')
    # start and end years to run for mean
    parser.add_argument(
        '--mean',
        '-m',
        metavar=('START', 'END'),
        type=int,
        nargs=2,
        default=[1980, 1995],
        help='Start and end year range for mean',
    )
    # fill value for ascii
    parser.add_argument(
        '--fill-value',
        '-f',
        type=float,
        default=np.nan,
        help='Output invalid value',
    )
    # print information about each input and output file
    parser.add_argument(
        '--verbose',
        '-V',
        default=False,
        action='store_true',
        help='Verbose output of run',
    )
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files',
    )
    # return the parser
    return parser


# Main program that calls gemb_smb_cumulative()
def main():
    parser = arguments()
    args = parser.parse_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logger = build_logger(__name__, level=loglevels[args.verbose])

    # run program
    gemb_smb_cumulative(
        args.infile, RANGE=args.mean, FILL_VALUE=args.fill_value, MODE=args.mode
    )


# run main program
if __name__ == '__main__':
    main()
