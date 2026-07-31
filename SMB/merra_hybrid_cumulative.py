#!/usr/bin/env python
"""
merra_hybrid_cumulative.py
Written by Tyler Sutterley (02/2023)
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
    Updated 02/2023: new doi for Medley (2022) Cryosphere paper
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

import re
import gzip
import uuid
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import model_harmonics as mdlhmc
from gravity_toolkit.utilities import build_logger

# dictionary defining the spatial reference system for each region
crs = dict(gris={}, ais={})
# Greenland (EPSG:3413)
crs['gris']['standard_name'] = 'Polar_Stereographic'
crs['gris']['grid_mapping_name'] = 'polar_stereographic'
crs['gris']['straight_vertical_longitude_from_pole'] = -45.0
crs['gris']['latitude_of_projection_origin'] = 90.0
crs['gris']['standard_parallel'] = 70.0
crs['gris']['scale_factor_at_projection_origin'] = 1.0
crs['gris']['false_easting'] = 0.0
crs['gris']['false_northing'] = 0.0
crs['gris']['semi_major_axis'] = 6378.137
crs['gris']['semi_minor_axis'] = 6356.752
crs['gris']['inverse_flattening'] = 298.257223563
crs['gris']['spatial_epsg'] = '3413'
# Antarctica (EPSG:3031)
crs['ais']['standard_name'] = 'Polar_Stereographic'
crs['ais']['grid_mapping_name'] = 'polar_stereographic'
crs['ais']['straight_vertical_longitude_from_pole'] = 0.0
crs['ais']['latitude_of_projection_origin'] = -90.0
crs['ais']['standard_parallel'] = -71.0
crs['ais']['scale_factor_at_projection_origin'] = 1.0
crs['ais']['false_easting'] = 0.0
crs['ais']['false_northing'] = 0.0
crs['ais']['semi_major_axis'] = 6378.137
crs['ais']['semi_minor_axis'] = 6356.752
crs['ais']['inverse_flattening'] = 298.257223563
crs['ais']['spatial_epsg'] = '3031'


# PURPOSE: calculate cumulative anomalies in MERRA-2 hybrid
# surface mass balance variables
def merra_hybrid_cumulative(
    base_dir, REGION, VERSION, RANGE=None, GZIP=False, MODE=0o775
):
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
    # get logger
    logger = logging.getLogger(__name__)

    # MERRA-2 hybrid directory
    DIRECTORY = base_dir.joinpath(VERSION)
    # set version parameters
    suffix = '.gz' if GZIP else ''
    if VERSION == 'v0':
        # input and output netCDF4 files
        args = (REGION.lower(), suffix)
        hybrid_file = 'm2_hybrid_p_minus_e_melt_{0}.nc{1}'.format(*args)
        output_file = 'm2_hybrid_cumul_{0}.nc{1}'.format(*args)
        # names of variables to read
        VARIABLES = ('p_minus_e', 'melt')
        AREA = None
        anomaly_flag = '_anomaly'
    elif VERSION in ('v1', 'v1.0'):
        # input and output netCDF4 files
        MAJOR_VERSION = re.match(r'((v\d+)(\.\d+)?)$', VERSION).group(2)
        args = (MAJOR_VERSION, REGION.lower(), suffix)
        hybrid_file = 'gsfc_fdm_smb_{0}_{1}.nc{2}'.format(*args)
        output_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
        # names of variables to read
        VARIABLES = ('runoff', 'rainfall', 'snowfall_minus_sublimation', 'SMB')
        AREA = None
        # flag to append to output netCDF4 variables
        anomaly_flag = '_anomaly'
    else:
        # input and output netCDF4 files
        FILE_VERSION = VERSION.replace('.', '_')
        args = (FILE_VERSION, REGION.lower(), suffix)
        firn_height_file = 'gsfc_fdm_{0}_{1}.nc{2}'.format(*args)
        hybrid_file = 'gsfc_fdm_smb_{0}_{1}.nc{2}'.format(*args)
        output_file = 'gsfc_fdm_smb_cumul_{0}_{1}.nc{2}'.format(*args)
        # names of variables to read
        VARIABLES = ('Me', 'Ra', 'Ru', 'Sn-Ev', 'SMB')
        AREA = 'iArea'
        # flag to append to output netCDF4 variables
        anomaly_flag = '_a'

    # Open the MERRA-2 Hybrid NetCDF file for reading
    input_file = DIRECTORY.joinpath(hybrid_file)
    if GZIP:
        # read as in-memory (diskless) netCDF4 dataset
        with gzip.open(str(input_file), mode='r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
    else:
        # read netCDF4 dataset
        fileID = netCDF4.Dataset(input_file, mode='r')

    # Output NetCDF file information
    logger.info(input_file)
    logger.info(list(fileID.variables.keys()))
    # dictionary describing the output netCDF4 structure
    struct = dict(
        dimensions=('time', 'x', 'y'),
        variables={},
    )
    # for each variable
    for v in VARIABLES:
        struct['variables'][f'{v}{anomaly_flag}'] = ('time', 'x', 'y')
    # output area variable if applicable
    if AREA:
        struct['variables'][AREA] = ('x', 'y')
    # dictionaries with data and attributes for output NetCDF file
    fd = {}
    attrs = dict(ROOT={})
    attributes_list = ['units', 'long_name', 'standard_name', 'comment']
    # global attributes of NetCDF file
    attrs['ROOT']['title'] = (
        f'Cumulative anomalies in GSFC-FDM{VERSION} variables '
        f'relative to {RANGE[0]:4d}-{RANGE[1]:4d}'
    )
    attrs['ROOT']['source'] = f'version {VERSION}'
    attrs['ROOT']['authors'] = 'Brooke Medley (NASA GSFC)'
    attrs['ROOT']['references'] = (
        'Medley, B., Neumann, T. A., Zwally, H. J., '
        'Smith, B. E., and Stevens, C. M.: Simulations of Firn Processes '
        'over the Greenland and Antarctic Ice Sheets: 1980--2021, '
        'The Cryosphere, https://doi.org/10.5194/tc-16-3971-2022, 2022.'
    )
    attrs['ROOT']['institution'] = 'NASA Goddard Space Flight Center (GSFC)'

    # list of input files for provenance
    lineage = []
    lineage.append(input_file.name)

    # create variable and attributes for projection
    fd['Polar_Stereographic'] = np.byte()
    # add projection attributes to dictionary
    attrs['Polar_Stereographic'] = crs[REGION]

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
    # Defining attributes for date
    attrs['time'] = dict(
        long_name='time, 5-daily resolution',
        units='decimal years, 5-daily resolution',
    )

    # extract areas from SMB or firn height file
    if AREA and (AREA in fileID.variables):
        fd[AREA] = fileID.variables[AREA][:].copy()
        # get each attribute for area variable if applicable
        attrs[AREA] = ncdf_attributes(fileID.variables[AREA], attributes_list)
        attrs[AREA]['grid_mapping'] = 'Polar_Stereographic'
    elif AREA:
        # Open the MERRA-2 Hybrid firn height file for reading
        input_firn_file = DIRECTORY.joinpath(firn_height_file)
        logger.debug(input_firn_file)
        lineage.append(input_firn_file.name)
        if GZIP:
            # read as in-memory (diskless) netCDF4 dataset
            with gzip.open(str(input_firn_file), mode='r') as f:
                fid1 = netCDF4.Dataset(uuid.uuid4().hex, memory=f.read())
        else:
            # read netCDF4 dataset
            fid1 = netCDF4.Dataset(input_file, mode='r')
        # copy area from firn height file
        fd[AREA] = fid1.variables[AREA][:].copy()
        # get each attribute for area variable if applicable
        attrs[AREA] = ncdf_attributes(fid1.variables[AREA], attributes_list)
        attrs[AREA]['grid_mapping'] = 'Polar_Stereographic'
        # close the firn height file
        fid1.close()

    # extract x and y coordinate arrays from grids if applicable
    # else create meshgrids of coordinate arrays
    if np.ndim(fileID.variables['x'][:]) == 2:
        xg = fileID.variables['x'][:].copy()
        yg = fileID.variables['y'][:].copy()
        fd['x'], fd['y'] = (xg[:, 0], yg[0, :])
    else:
        fd['x'] = fileID.variables['x'][:].copy()
        fd['y'] = fileID.variables['y'][:].copy()
        xg, yg = np.meshgrid(fd['x'], fd['y'], indexing='ij')
    # time is year decimal at time step 5 days
    time_step = 5.0 / 365.25
    # calculate mean period for MERRA-2
    (tt,) = np.nonzero((fd['time'] >= RANGE[0]) & (fd['time'] < (RANGE[1] + 1)))

    # for each variable
    for v in VARIABLES:
        # append anomaly flag
        var = f'{v}{anomaly_flag}'
        # copy data and remove singleton dimensions
        # set masks
        DATA = np.ma.masked_equal(
            fileID.variables[v][:].squeeze(),
            fileID.variables[v]._FillValue,
        )
        # get each attribute for variable if applicable
        attrs[var] = ncdf_attributes(
            fileID.variables[v], attributes_list, replace=[(' per year', '')]
        )
        attrs[var]['grid_mapping'] = 'Polar_Stereographic'
        # input shape of MERRA-2 Hybrid firn data
        nt, nx, ny = np.shape(DATA)

        # cumulative mass anomalies calculated by removing mean balance flux
        # mean of data for variable (converted from yearly rate)
        MEAN = np.mean(DATA.data[tt, :, :] * time_step, axis=0)
        # indices of specified ice mask at the first slice
        i, j = np.nonzero(~DATA.mask[0, :, :])
        valid_count = np.count_nonzero(~DATA.mask[0, :, :])
        # allocate for output variable
        fd[var] = np.ma.zeros((nt, nx, ny), fill_value=DATA.fill_value)
        fd[var].mask = DATA.mask | np.isnan(DATA.data)
        CUMULATIVE = np.zeros((valid_count))
        # calculate output cumulative anomalies for variable
        for t in range(nt):
            # convert mass flux from yearly rate and
            # calculate cumulative anomalies at time t
            CUMULATIVE += DATA.data[t, i, j] * time_step - MEAN[i, j]
            fd[var].data[t, i, j] = CUMULATIVE.copy()
        # replace masked values with fill value
        fd[var].data[fd[var].mask] = fd[var].fill_value
    # close the NetCDF files
    fileID.close()

    # Output NetCDF filename
    output_cumulative_file = DIRECTORY.joinpath(output_file)
    logger.info(output_cumulative_file)

    # output MERRA-2 data file with cumulative data
    if GZIP:
        # write data to in-memory netCDF4 file
        nc_buffer = mdlhmc.spatial.to_netCDF4(
            uuid.uuid4().hex,
            fd,
            attrs,
            struct,
            mode='w',
            memory=True,
        )
        # write MERRA-2 data file to gzipped file
        with gzip.open(str(output_file), 'wb') as f:
            f.write(nc_buffer)
    else:
        # write data to netCDF4 file
        mdlhmc.spatial.to_netCDF4(
            output_file,
            fd,
            attrs,
            struct,
            mode='w',
        )
    # change the permissions mode
    output_file.chmod(mode=MODE)
    # change the permissions mode
    output_cumulative_file.chmod(mode=MODE)


# PURPOSE: get attributes for netCDF4 files and variable
def ncdf_attributes(nc, attributes_list, replace=[]):
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
        # replace attribute values if applicable
        for key, val in replace:
            att_val = att_val.replace(key, val)
        # add attribute to dictionary
        attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes


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
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # region of firn model
    parser.add_argument(
        '--region',
        '-R',
        type=str,
        default='gris',
        choices=['gris', 'ais'],
        help='Region of firn model to calculate',
    )
    # version of firn model
    versions = ['v0', 'v1', 'v1.0', 'v1.1', 'v1.2', 'v1.2.1']
    parser.add_argument(
        '--version',
        '-v',
        type=str,
        default='v1.2.1',
        choices=versions,
        help='Version of firn model to calculate',
    )
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
    # netCDF4 files are gzip compressed
    parser.add_argument(
        '--gzip',
        '-G',
        default=False,
        action='store_true',
        help='netCDF4 file is locally gzip compressed',
    )
    # print information about each input and output file
    parser.add_argument(
        '--verbose',
        '-V',
        action='count',
        default=0,
        help='Verbose output of processing run',
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


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # create logger
    loglevels = [logging.CRITICAL, logging.INFO, logging.DEBUG]
    logger = build_logger(__name__, level=loglevels[args.verbose])

    # run program
    merra_hybrid_cumulative(
        args.directory,
        args.region,
        args.version,
        RANGE=args.mean,
        GZIP=args.gzip,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
