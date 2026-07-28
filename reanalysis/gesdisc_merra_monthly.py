#!/usr/bin/env python
"""
gesdisc_merra_monthly.py
Written by Tyler Sutterley (07/2026)

Creates monthly MERRA-2 3D model level products syncing data from the
    Goddard Earth Sciences Data and Information Server Center (GES DISC)
    https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
Combines daily model level outputs into monthly averages

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_merra_monthly.py --user <username>
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: Working data directory
    -s X, --shortname X: MERRA-2 product shortname
    -v X, --version X: MERRA-2 version
    -Y X, --year X: years to sync
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -V, --verbose: Output information for each output file
    -M X, --mode X: Local permissions mode of the files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 07/2026: use numpy summation to calculate daily averages
        added axis to variable attributes to attempt to retrieve
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 03/2023: use full path to output file in verbose logging
    Updated 12/2022: single implicit import of spherical harmonic tools
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: add hard coded GES DISC subsetting api host option
    Updated 08/2022: adjust time range for CMR queries
    Updated 06/2022: use CMR queries to find reanalysis granules
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added option for connection timeout (in seconds)
        use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 12/2020: use argparse to set command line parameters
        using utilities program to build opener
    Updated 09/2019: added ssl context to urlopen headers
        modified regular expression to check if suffix of file is nc4
    Updated 08/2019: new GESDISC server and links list file format
    Updated 06/2018: using python3 compatible octal, input and urllib
    Written 03/2018
"""

from __future__ import print_function

import sys
import os
import re
import time
import netrc
import getpass
import logging
import netCDF4
import pathlib
import argparse
import builtins
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc


# PURPOSE: sync local MERRA-2 files with GESDISC server
def gesdisc_merra_monthly(
    base_dir,
    SHORTNAME,
    VERSION=None,
    YEARS=None,
    TIMEOUT=None,
    LOG=False,
    VERBOSE=False,
    MODE=None,
):
    # full path to MERRA-2 directory
    base_dir = pathlib.Path(base_dir).expanduser().absolute()
    # check if DIRECTORY exists and recursively create if not
    DIRECTORY = base_dir.joinpath('MERRA-2')
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    if LOG:
        # format: NASA_GESDISC_MERRA2_monthly_2002-04-01.log
        today = time.strftime('%Y-%m-%d', time.localtime())
        output_logfile = f'NASA_GESDISC_MERRA2_monthly_{today}.log'
        LOGFILE = DIRECTORY.joinpath(output_logfile)
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logging.basicConfig(stream=fid, level=loglevel)
        logging.info(f'NASA MERRA-2 Sync Log ({today})')
        PRODUCT = f'{SHORTNAME}.{VERSION}'
        logging.info(f'PRODUCT: {PRODUCT}')
    else:
        # standard output (terminal output)
        fid = sys.stdout
        logging.basicConfig(stream=fid, level=loglevel)

    # regular expression for grouping months from daily data
    regex_pattern = r'MERRA2_(\d+).(.*?).(\d{4})(\d{2})(\d{2})(.*?).nc[4]?$'
    rx1 = re.compile(regex_pattern, re.VERBOSE)

    # output variable names
    VARNAME = 'PS'
    TNAME = 'T'
    QNAME = 'QV'
    LONNAME = 'lon'
    LATNAME = 'lat'
    LEVELNAME = 'lev'
    TIMENAME = 'time'
    # output dimensions
    ntime, nlevels, nlat, nlon = (1, 72, 361, 576)
    # dictionary defining output structure
    struct = dict(
        dimensions=(TIMENAME, LEVELNAME, LATNAME, LONNAME),
        variables={
            VARNAME: (TIMENAME, LATNAME, LONNAME),
            TNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
            QNAME: (TIMENAME, LEVELNAME, LATNAME, LONNAME),
        },
        shape={
            VARNAME: (ntime, nlat, nlon),
            TNAME: (ntime, nlevels, nlat, nlon),
            QNAME: (ntime, nlevels, nlat, nlon),
        },
    )
    # dictionary defining file-level and variable attributes
    attributes = dict(ROOT={})
    # file-level attributes to retrieve
    root_attributes = [
        'Contact',
        'Conventions',
        'Institution',
        'References',
        'Format',
        'SpatialCoverage',
        'VersionID',
        'identifier_product_doi_authority',
        'identifier_product_doi',
        'ShortName',
        'LongName',
        'Title',
        'DataResolution',
        'LatitudeResolution',
        'LongitudeResolution',
        'SouthernmostLatitude',
        'NorthernmostLatitude',
        'WesternmostLongitude',
        'EasternmostLongitude',
    ]
    # variable attributes to retrieve
    variable_attributes = [
        'axis',
        'calendar',
        'long_name',
        'positive',
        'standard_name',
        'units',
        'valid_range',
    ]

    # for each unique date
    for YEAR in YEARS:
        dpm = gravtk.time.calendar_days(YEAR)
        # for each month of the year
        for i, days_per_month in enumerate(dpm):
            # year and month as strings
            YY = f'{YEAR:4d}'
            MM = f'{i + 1:02d}'
            # start and end date for query
            start_date = f'{YY}-{MM}-{1:02.0f}T00:00:00'
            end_date = f'{YY}-{MM}-{days_per_month:02.0f}T23:59:59'
            # output time units
            time_units = f'minutes since {YY}-{MM}-01 00:00:00'
            # query for data
            ids, urls, mtimes = mdlhmc.utilities.cmr(
                SHORTNAME,
                version=VERSION,
                start_date=start_date,
                end_date=end_date,
                provider='GES_DISC',
                verbose=VERBOSE,
            )
            # skip years and months without any data
            if not ids:
                continue
            # dictionary with output data
            dinput = {}
            # dictionary with count for converting totals to means
            count = {}
            # allocate arrays for each variable and count of valid values
            dinput[TIMENAME] = np.zeros((ntime))
            count[TIMENAME] = np.zeros((ntime), dtype='i')
            for var, shape in struct['shape'].items():
                dinput[var] = np.zeros(shape)
                count[var] = np.zeros(shape, dtype='i')
            # for each url
            for id, url, mtime in zip(ids, urls, mtimes):
                # build subsetting API url for granule
                request_url = mdlhmc.utilities.build_request(
                    SHORTNAME,
                    VERSION,
                    url,
                    host='https://goldsmr5.gesdisc.eosdis.nasa.gov',
                    variables=struct['variables'].keys(),
                )
                # Create and submit request. There are a wide range of exceptions
                # that can be thrown here, including HTTPError and URLError.
                response = mdlhmc.utilities.from_http(
                    request_url,
                    timeout=TIMEOUT,
                    context=None,
                    verbose=VERBOSE,
                    fid=fid,
                )
                response.seek(0)
                # open remote file with netCDF4
                fileID = netCDF4.Dataset(id, 'r', memory=response.read())
                MOD, DATASET, Y, M, D, AUX = rx1.findall(id).pop()
                # extract dimension variables
                for dim in (LEVELNAME, LATNAME, LONNAME):
                    dinput[dim] = fileID.variables[dim][:].copy()
                    # extract variable attributes
                    attributes[dim] = ncdf_attributes(
                        fileID[dim], variable_attributes
                    )
                # extract time dimension variable
                TIME = fileID.variables[TIMENAME]
                # get the epoch and units from the attributes
                epoch1, to_secs = gravtk.time.parse_date_string(TIME.units)
                epoch2, _ = gravtk.time.parse_date_string(time_units)
                # convert delta times to new epoch
                delta_time = gravtk.time.convert_delta_time(
                    to_secs * TIME[:],
                    epoch1=epoch1,
                    epoch2=epoch2,
                    scale=1.0 / to_secs,
                )
                # add over time slices products to monthly output
                dinput[TIMENAME][0] += np.sum(delta_time)
                count[TIMENAME][0] += len(delta_time)
                # extract variable attributes
                attributes[TIMENAME] = ncdf_attributes(
                    fileID[TIMENAME], variable_attributes
                )
                # update time units attribute
                attributes['time']['units'] = time_units
                # for each variable
                for var in struct['variables'].keys():
                    # surface pressure
                    # air temperature
                    # specific humidity
                    tmp = fileID.variables[var][:]
                    valid = tmp != fileID.variables[var]._FillValue
                    tmp[~valid] = 0.0
                    # calculate sum over time axis
                    dinput[var][0, ...] += tmp.sum(axis=0)
                    count[var][0, ...] += valid.sum(axis=0)
                    # extract variable attributes
                    attributes[var] = ncdf_attributes(
                        fileID[var], variable_attributes
                    )
                # get the fill value for the surface variables
                fill_value = fileID.variables[VARNAME]._FillValue
                # try to get root attributes
                attributes['ROOT'] = ncdf_attributes(fileID, root_attributes)
                # close the input file from remote url
                fileID.close()
            # calculate mean from totals
            dinput[TIMENAME] /= count[TIMENAME]
            for var in struct['variables'].keys():
                # find valid values
                valid = count[var] > 0
                valid_indices = np.nonzero(valid)
                # calculate mean fields from totals
                dinput[var][valid_indices] /= count[var][valid_indices]
                # replace points where no values with fill_value
                complementary_indices = np.nonzero(~valid)
                dinput[var][complementary_indices] = fill_value
            # output to netCDF4 file (replace hour variable with monthly)
            DATASET = update_attribute(DATASET)
            FILE = f'MERRA2_{MOD}.{DATASET}.{YY}{MM}.SUB.nc'
            local_file = DIRECTORY.joinpath(FILE)
            ncdf_model_write(
                dinput,
                attributes,
                fill_value,
                struct,
                FILENAME=local_file,
            )
            # keep remote modification time of file and local access time
            os.utime(local_file, (local_file.stat().st_atime, mtime))
            # set permissions mode to MODE
            local_file.chmod(mode=MODE)

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
            attributes[att_name] = update_attribute(att_val)
    # return the dictionary of attributes
    return attributes


# PURPOSE: change attributes from inst3_3d to tavgM_3d
def update_attribute(name):
    output = re.sub(r'inst\d+_3d', r'tavgM_3d', name, re.I)
    return output


# PURPOSE: write output model layer fields data to file
def ncdf_model_write(
    output,
    attributes,
    fill_value,
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
        nc[var] = fileID.createVariable(
            var,
            output[var].dtype,
            dimensions,
            fill_value=fill_value,
            zlib=True,
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
        description="""Creates monthly MERRA-2 3D model level
            products syncing data from the Goddard Earth Sciences
            Data and Information Server Center (GES DISC)
            """
    )
    # command line parameters
    # NASA Earthdata credentials
    parser.add_argument(
        '--user',
        '-U',
        type=str,
        default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login',
    )
    parser.add_argument(
        '--password',
        '-W',
        type=str,
        default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login',
    )
    parser.add_argument(
        '--netrc',
        '-N',
        type=pathlib.Path,
        default=pathlib.Path.home().joinpath('.netrc'),
        help='Path to .netrc file for authentication',
    )
    # working data directory
    parser.add_argument(
        '--directory',
        '-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory',
    )
    # MERRA-2 product shortname
    parser.add_argument(
        '--shortname',
        '-s',
        type=str,
        default='M2I3NVASM',
        help='MERRA-2 product shortname',
    )
    # MERRA-2 version
    parser.add_argument(
        '--version',
        '-v',
        type=str,
        default='5.12.4',
        help='MERRA-2 version',
    )
    # years to download
    now = time.gmtime()
    parser.add_argument(
        '--year',
        '-Y',
        type=int,
        nargs='+',
        default=range(2000, now.tm_year + 1),
        help='Years of model outputs to sync',
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
    # NASA_GESDISC_MERRA2_monthly_2002-04-01.log
    parser.add_argument(
        '--log',
        '-l',
        default=False,
        action='store_true',
        help='Output log file',
    )
    # print information about each output file
    parser.add_argument(
        '--verbose',
        '-V',
        default=False,
        action='store_true',
        help='Verbose output of run',
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

    # NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    # get NASA Earthdata credentials
    try:
        args.user, _, args.password = netrc.netrc(args.netrc).authenticators(
            URS
        )
    except:
        # check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = f'Username for {URS}: '
            args.user = builtins.input(prompt)
        # enter password securely from command-line
        if not args.password:
            prompt = f'Password for {args.user}@{URS}: '
            args.password = getpass.getpass(prompt)

    # build a urllib opener for NASA GESDISC
    # Add the username and password for NASA Earthdata Login system
    mdlhmc.utilities.build_opener(
        args.user,
        args.password,
        password_manager=True,
        authorization_header=False,
    )

    # check internet connection before attempting to run program
    HOST = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/'
    if mdlhmc.utilities.check_credentials(HOST):
        gesdisc_merra_monthly(
            args.directory,
            args.shortname,
            VERSION=args.version,
            YEARS=args.year,
            TIMEOUT=args.timeout,
            LOG=args.log,
            VERBOSE=args.verbose,
            MODE=args.mode,
        )


# run main program
if __name__ == '__main__':
    main()
