#!/usr/bin/env python
"""
ucar_gdex_jra3q_levels.py
Written by Tyler Sutterley (07/2026)

Downloads JRA-3Q model-level analysis products provided by the
    NCAR/UCAR Geoscience Data Exchange (GDEX)

JRA-3Q: Japanese Reanalysis for Three Quarters of a Century
    https://gdex.ucar.edu/datasets/d640000

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: Working data directory
    --variable X, -V X: Variable to download
        * hgt-hyb: geopotential height
    -Y X, --year X: Years to download
    -I, --invariant: Retrieve the model invariant parameters
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: Output log of files downloaded
    -M X, --mode=X: Permission mode of directories and files downloaded

PROGRAM DEPENDENCIES:
    spatial.py: spatial data class for reading, writing and processing data
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Written 07/2026
"""

from __future__ import print_function

import sys
import os
import io
import re
import copy
import time
import logging
import netCDF4
import pathlib
import argparse
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local JRA-3Q files with UCAR/NCAR GDEX server
def ucar_gdex_download(
    base_dir,
    VARIABLE,
    YEARS=None,
    TIMEOUT=None,
    LOG=False,
    MODE=None,
):
    # directory setup
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    DIRECTORY = base_dir.joinpath('JRA-3Q')
    # check if log directory exists and recursively create if not
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    if LOG:
        # format: UCAR_GDEX_JRA-3Q_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'UCAR_GDEX_JRA-3Q_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=logging.INFO)
        logging.info(f'UCAR JRA-3Q Sync Log ({today})')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=logging.INFO)

    # UCAR GDEX host
    HOST = 'https://gdex.ucar.edu/'
    dataset_id = 'd640000'
    product_id = '3'
    PRODUCT = 'anl_mdl'
    TIMENAME = 'time'
    LEVELNAME = 'hybrid_level'
    HALFNAME = 'hybrid_half_level'
    LATNAME = 'lat'
    LONNAME = 'lon'
    # netCDF4 variable name for the requested product
    VARNAME = f'{VARIABLE}-an-gauss'
    POINTS = 'original_number_of_grid_points_per_latitude_circle'
    WEIGHT = 'weight'
    ANAME, BNAME = ('a_hybrid_half_level', 'b_hybrid_half_level')
    AINTERFACE, BINTERFACE = ('a_hybrid_level', 'b_hybrid_level')
    # build a file filter for the product and variable
    query_params = {}
    query_params['filter_wfile'] = VARIABLE
    kwargs = '?' + gravtk.utilities.urlencode(query_params)

    # dimension sizes for the output netCDF4 file
    ntime, nlevels, nlat, nlon = (1, 100, 480, 960)
    # dictionary defining output structure
    struct = dict(
        dimensions=(TIMENAME, LEVELNAME, HALFNAME, LATNAME, LONNAME),
        variables={
            VARNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
            POINTS: (LATNAME,),
            WEIGHT: (LATNAME,),
            ANAME: (HALFNAME,),
            BNAME: (HALFNAME,),
            AINTERFACE: (LEVELNAME,),
            BINTERFACE: (LEVELNAME,),
        },
        shape={
            VARNAME: (ntime, nlevels, nlat, nlon),
            POINTS: (nlat,),
            WEIGHT: (nlat,),
            ANAME: (nlevels + 1,),
            BNAME: (nlevels + 1,),
            AINTERFACE: (nlevels,),
            BINTERFACE: (nlevels,),
        },
    )
    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    # file-level attributes to retrieve
    root_attributes = [
        'jma_data_provider',
        'jma_data_provider_url',
        'jma_data_original_grid',
        'jma_data_title',
        'jma_data_url',
        'jma_primary_publication',
        'jma_primary_publication_url',
        'rda_dataset',
        'rda_dataset_title',
        'rda_dataset_url',
        'rda_dataset_doi',
        'rda_dataset_license',
        'rda_data_output_grid',
        'Conventions',
    ]
    # variable attributes to retrieve
    variable_attributes = [
        'axis',
        'calendar',
        'data_type',
        'formula',
        'long_name',
        'jma_short_name',
        'short_name',
        'standard_name',
        'time_step',
        'units',
        'valid_range',
        'shape_of_the_earth_table'
        'shape_of_the_earth_number'
        'shape_of_the_earth_meaning'
        'stephan_boltzmann_constant'
        'earth_radius'
        'angular_speed_of_earth_rotation'
        'gravitational_accleration'
        'gas_constant_for_dry_air'
        'specific_heat_of_dry_air_at_constant_pressure'
        'latent_heat_of_vaporization'
        'solar_constant',
    ]

    # find data directories for year
    dirs, _ = mdlhmc.utilities.ucar_list(
        [HOST, 'datasets', dataset_id, 'filelist', product_id],
        tdclass='Group Name',
        timeout=TIMEOUT,
    )
    # model years in each directory
    years = [re.findall(rf'(\d+){product_id.zfill(2)}', d).pop() for d in dirs]
    # reduce list of directories to only those for the requested years
    if YEARS is not None:
        directories = copy.copy(dirs)
        dirs = [directories[years.index(y)] for y in YEARS if y in years]

    # for each directory in the (reduced) list of directories
    for d in dirs:
        (year,) = re.findall(rf'(\d+){product_id.zfill(2)}', d)
        # for each month
        for m in range(12):
            # find data files for year and month
            dirparts = mdlhmc.utilities.url_split(d)
            cols, mods = mdlhmc.utilities.ucar_list(
                [HOST, *dirparts, kwargs],
                timeout=TIMEOUT,
                pattern=f'{year}{str(m + 1).zfill(2)}',
                sort=True,
            )
            # skip to next month if no files found
            if not any(cols):
                continue
            # dictionary with output data
            dinput = {}
            # dictionary with count for converting totals to means
            count = {}
            # allocate arrays for each variable and count of valid values
            dinput[TIMENAME] = np.zeros((ntime))
            count[TIMENAME] = np.zeros((ntime), dtype='i')
            shape = struct['shape'][VARNAME]
            dinput[VARNAME] = np.zeros(shape)
            count[VARNAME] = np.zeros(shape, dtype='i')
            # for each file in the list of files
            for colname, collastmod in zip(cols, mods):
                # download file from UCAR GDEX server
                logging.debug(colname)
                fileparts = mdlhmc.utilities.url_split(colname)
                # create and submit request for file
                response = mdlhmc.utilities.from_http(
                    colname,
                    timeout=TIMEOUT,
                    local=None,
                    verbose=True,
                    fid=fid,
                    mode=MODE,
                )
                response.seek(0)
                # open remote file with netCDF4
                fileID = netCDF4.Dataset(id, 'r', memory=response.read())
                # extract dimension variables
                for dim in (LEVELNAME, HALFNAME, LATNAME, LONNAME):
                    dinput[dim] = fileID.variables[dim][:].copy()
                    # extract variable attributes
                    attributes[dim] = ncdf_attributes(
                        fileID[dim], variable_attributes
                    )
                # extract gaussian grid variables
                for var in (POINTS, WEIGHT):
                    dinput[var] = fileID.variables[var][:].copy()
                    # extract variable attributes
                    attributes[var] = ncdf_attributes(
                        fileID[var], variable_attributes
                    )
                # extract level variables
                for var in (ANAME, BNAME, AINTERFACE, BINTERFACE):
                    dinput[var] = fileID.variables[var][:].copy()
                    # extract variable attributes
                    attributes[var] = ncdf_attributes(
                        fileID[var], variable_attributes
                    )
                # extract time dimension variable
                delta_time = fileID.variables[TIMENAME][:].copy()
                # add over time slices products to monthly output
                dinput[TIMENAME][0] += np.sum(delta_time)
                count[TIMENAME][0] += len(delta_time)
                # extract variable attributes
                attributes[TIMENAME] = ncdf_attributes(
                    fileID[TIMENAME], variable_attributes
                )
                # variables of interest for the requested product
                for var in (VARNAME,):
                    # geopotential height
                    tmp = fileID.variables[var][:]
                    # replace invalid values with 0.0
                    valid = tmp != fileID[var]._FillValue
                    tmp[~valid] = 0.0
                    # calculate sum over time axis
                    dinput[var][0, ...] += tmp.sum(axis=0)
                    count[var][0, ...] += valid.sum(axis=0)
                    # extract variable attributes
                    attributes[var] = ncdf_attributes(
                        fileID[var], variable_attributes
                    )
                # get the fill value for the output variables
                fill_value = fileID[VARNAME]._FillValue
                # try to get root attributes
                attributes['ROOT'] = ncdf_attributes(fileID, root_attributes)
                # close the input file from remote url
                fileID.close()

            # calculate mean time
            dinput[TIMENAME] /= count[TIMENAME]
            # calculate mean of each variable
            for var in (VARNAME,):
                # find valid values
                count[var] = np.ma.masked_equal(count[var], 0)
                # calculate mean fields from totals
                dinput[var] = dinput[var] / count[var]
                dinput[var].set_fill_value(fill_value)
                dinput[var].data[dinput[var].mask] = fill_value

            # convert time to strings
            (datetime,) = gravtk.time.to_string(
                dinput[TIMENAME],
                attributes[TIMENAME]['units'],
                strftime='%Y%m',
            )

            # output filename
            filename = f'jra3q.{PRODUCT}.{VARNAME}.{datetime}.nc'
            output = DIRECTORY.joinpath(filename)
            ncdf_model_write(
                dinput,
                attributes,
                struct,
                FILENAME=output,
            )
            # keep remote modification time of file and local access time
            os.utime(output, (output.stat().st_atime, collastmod))
            # set permissions mode to MODE
            output.chmod(mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        fid.close()
        LOGFILE.chmod(mode=MODE)


# PURPOSE: get attributes for a variable
def ncdf_attributes(nc, attributes_list):
    # output dictionary of attributes for variable
    attributes = {}
    # for each attribute to try to get
    for att_name in attributes_list:
        try:
            att_val = nc.getncattr(att_name)
        except Exception as exc:
            logging.debug(f'Attribute {att_name} not found in {nc.name}')
            pass
        else:
            attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes


# PURPOSE: write output model layer fields data to file
def ncdf_model_write(
    output,
    attributes,
    struct,
    FILENAME=None,
):
    # opening NetCDF4 file for writing
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    fileID = netCDF4.Dataset(FILENAME, 'w', format='NETCDF4')

    # dictionary with NetCDF4 variable objects
    nc = {}
    # defining the NetCDF4 dimensions
    for dim in struct['dimensions']:
        fileID.createDimension(dim, len(output[dim]))
        nc[dim] = fileID.createVariable(dim, output[dim].dtype, (dim,))
        # add data to NetCDF4 dimension variable
        nc[dim][:] = output[dim].copy()
        # set netCDF4 attributes for dimensions
        for att_name, att_val in attributes[dim].items():
            nc[dim].setncattr(att_name, att_val)

    # defining the NetCDF4 variables
    for var, dimensions in struct['variables'].items():
        if hasattr(output[var], 'fill_value'):
            nc[var] = fileID.createVariable(
                var,
                output[var].dtype,
                dimensions,
                fill_value=output[var].fill_value,
                zlib=True,
            )
        else:
            nc[var] = fileID.createVariable(
                var,
                output[var].dtype,
                dimensions,
            )
        # add data to NetCDF4 variable
        nc[var][:] = output[var].copy()
        # set netCDF4 attributes for variables
        for att_name, att_val in attributes[var].items():
            nc[var].setncattr(att_name, att_val)

    # Defining file-level attributes
    for att_name, att_val in attributes['ROOT'].items():
        fileID.setncattr(att_name, att_val)
    # add attribute for date created
    fileID.date_created = time.strftime('%Y-%m-%d', time.localtime())
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    fileID.reference = f'Output from {pathlib.Path(sys.argv[0]).name}'

    # Output NetCDF structure information
    logging.info(str(FILENAME))
    logging.info(list(fileID.variables.keys()))

    # Closing the NetCDF file
    fileID.close()


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Downloads JRA-3Q model-level analysis products
            provided by the NCAR/UCAR Geoscience Data Exchange (GDEX)
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
    # JRA-3Q model-level analysis variable
    parser.add_argument(
        '--variable',
        '-v',
        type=str,
        default='hgt-hyb',
        help='JRA-3Q model-level analysis variable',
    )
    # model years to download
    parser.add_argument(
        '--year',
        '-Y',
        type=str,
        nargs='+',
        help='Years to download',
    )
    # connection timeout
    parser.add_argument(
        '--timeout',
        '-t',
        type=int,
        default=360,
        help='Timeout in seconds for blocking operations',
    )
    # Output log file in form
    # UCAR_GDEX_JRA-3Q_2002-04-01.log
    parser.add_argument(
        '--log',
        '-l',
        default=False,
        action='store_true',
        help='Output log file',
    )
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument(
        '--mode',
        '-M',
        type=lambda x: int(x, base=8),
        default=0o775,
        help='Permission mode of directories and files synced',
    )
    # return the parser
    return parser


# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args, _ = parser.parse_known_args()

    # download JRA-3Q products
    ucar_gdex_download(
        args.directory,
        args.variable,
        YEARS=args.year,
        TIMEOUT=args.timeout,
        LOG=args.log,
        MODE=args.mode,
    )


# run main program
if __name__ == '__main__':
    main()
