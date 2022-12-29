#!/usr/bin/env python
u"""
gesdisc_merra_monthly.py
Written by Tyler Sutterley (12/2022)

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
import argparse
import builtins
import numpy as np
import gravity_toolkit as gravtk
import model_harmonics as mdlhmc

# PURPOSE: sync local MERRA-2 files with GESDISC server
def gesdisc_merra_monthly(base_dir, SHORTNAME, VERSION=None, YEARS=None,
    TIMEOUT=None, LOG=False, VERBOSE=False, MODE=None):
    # full path to MERRA-2 directory
    DIRECTORY = os.path.join(base_dir,'MERRA-2')
    # check if DIRECTORY exists and recursively create if not
    if not os.access(os.path.join(DIRECTORY), os.F_OK):
        os.makedirs(os.path.join(DIRECTORY), mode=MODE, exist_ok=True)

    # create log file with list of synchronized files (or print to terminal)
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    if LOG:
        # format: NASA_GESDISC_MERRA2_monthly_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = f'NASA_GESDISC_MERRA2_monthly_{today}.log'
        fid = open(os.path.join(DIRECTORY,LOGFILE), mode='w', encoding='utf8')
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
    # time format for CMR queries
    isotime_format = '{0}-{1}-{2:02.0f}T{3:02.0f}:{4:02.0f}:{5:02.0f}'
    # output dimensions
    nlevels,nlat,nlon = (72,361,576)
    # dictionary of variable attributes
    attributes = dict(ROOT={})
    # file-level attributes to retrieve
    root_attributes = ['Contact', 'Conventions', 'Institution',
        'References', 'Format', 'SpatialCoverage', 'VersionID',
        'identifier_product_doi_authority', 'identifier_product_doi',
        'ShortName', 'LongName', 'Title', 'DataResolution',
        'LatitudeResolution', 'LongitudeResolution',
        'SouthernmostLatitude', 'NorthernmostLatitude',
        'WesternmostLongitude', 'EasternmostLongitude']

    # for each unique date
    for YEAR in YEARS:
        dpm = gravtk.time.calendar_days(YEAR)
        # for each month of the year
        for i,days_per_month in enumerate(dpm):
            # year and month as strings
            YY = f'{YEAR:4d}'
            MM = f'{i+1:02d}'
            # start and end date for query
            start_date = isotime_format.format(YY,MM,1.0,0.0,0.0,0.0)
            end_date = isotime_format.format(YY,MM,days_per_month,23.0,59.0,59.0)
            # query for data
            ids,urls,mtimes = mdlhmc.utilities.cmr(SHORTNAME, version=VERSION,
                start_date=start_date, end_date=end_date,
                provider='GES_DISC', verbose=VERBOSE)
            # skip years and months without any data
            if not ids:
                continue
            # python dictionary with output data
            dinput = {}
            dinput[TIMENAME] = np.zeros((1))
            dinput[VARNAME] = np.zeros((1,nlat,nlon))
            dinput[TNAME] = np.zeros((1,nlevels,nlat,nlon))
            dinput[QNAME] = np.zeros((1,nlevels,nlat,nlon))
            # python dictionary with count for converting totals to means
            count = {}
            count[TIMENAME] = np.zeros((1))
            count[VARNAME] = np.zeros((1,nlat,nlon))
            count[TNAME] = np.zeros((1,nlevels,nlat,nlon))
            count[QNAME] = np.zeros((1,nlevels,nlat,nlon))
            # for each url
            for id,url,mtime in zip(ids,urls,mtimes):
                # build subsetting API url for granule
                request_url = mdlhmc.utilities.build_request(SHORTNAME,
                    VERSION, url, host='https://goldsmr5.gesdisc.eosdis.nasa.gov',
                    variables=[VARNAME,TNAME,QNAME])
                # Create and submit request. There are a wide range of exceptions
                # that can be thrown here, including HTTPError and URLError.
                response = mdlhmc.utilities.from_http(request_url,
                    timeout=TIMEOUT,
                    context=None,
                    verbose=VERBOSE,
                    fid=fid)
                response.seek(0)
                # open remote file with netCDF4
                fileID = netCDF4.Dataset(id,'r',memory=response.read())
                MOD,DATASET,Y,M,D,AUX = rx1.findall(id).pop()
                # extract dimension variables
                nt, = fileID.variables[TIMENAME].shape
                for dim in (LEVELNAME,LATNAME,LONNAME):
                    dinput[dim] = fileID.variables[dim][:].copy()
                    # extract variable attributes
                    attributes[dim] = ncdf_attributes(fileID,dim)
                # bad value
                fill_value = fileID.variables[VARNAME]._FillValue
                # add over time slices products to monthly output
                for t in range(nt):
                    TIME = fileID.variables[TIMENAME][t].astype('f')
                    dinput[TIMENAME][0] += TIME
                    count[TIMENAME][0] += 1.0
                    # surface pressure
                    PS = fileID.variables[VARNAME][t,:,:].copy()
                    ii,jj = np.nonzero(PS != fill_value)
                    dinput[VARNAME][0,ii,jj] += PS[ii,jj]
                    count[VARNAME][0,ii,jj] += 1.0
                    # air temperature
                    T = fileID.variables[TNAME][t,:,:,:].copy()
                    ii,jj,kk = np.nonzero(T != fill_value)
                    dinput[TNAME][0,ii,jj,kk] += T[ii,jj,kk]
                    count[TNAME][0,ii,jj,kk] += 1.0
                    # specific humidity
                    QV = fileID.variables[QNAME][t,:,:,:].copy()
                    ii,jj,kk = np.nonzero(QV != fill_value)
                    dinput[QNAME][0,ii,jj,kk] += QV[ii,jj,kk]
                    count[QNAME][0,ii,jj,kk] += 1.0
                    # get attributes for each variable
                    for var in (TIMENAME, VARNAME, TNAME, QNAME):
                        # extract variable attributes
                        attributes[var] = ncdf_attributes(fileID,var)
                # get each root attribute of interest
                for att_name in root_attributes:
                    try:
                        att_val = fileID.getncattr(att_name)
                        att_val = re.sub(r'inst\d+_3d',r'instM_3d',att_val)
                    except Exception as e:
                        pass
                    else:
                        attributes['ROOT'][att_name] = att_val
                # close the input file from remote url
                fileID.close()
            # calculate mean from totals
            dinput[TIMENAME] /= count[TIMENAME]
            for key in [VARNAME,TNAME,QNAME]:
                # find valid values
                valid_indices = np.nonzero(count[key] > 0)
                dinput[key][valid_indices] /= count[key][valid_indices]
                # replace points where no values with fill_value
                complementary_indices = np.nonzero(count[key] == 0)
                dinput[key][complementary_indices] = fill_value
            # output to netCDF4 file (replace hour variable with monthly)
            DATASET = re.sub(r'inst\d+_3d',r'instM_3d',DATASET)
            attributes['time']['units'] = f'minutes since {YY}-{MM}-01 00:00:00'
            local_file = f'MERRA2_{MOD}.{DATASET}.{YY}{MM}.SUB.nc'
            ncdf_model_write(dinput, attributes, fill_value,
                VARNAME=VARNAME, TNAME=TNAME, QNAME=QNAME,
                LONNAME=LONNAME, LATNAME=LATNAME, LEVELNAME=LEVELNAME,
                TIMENAME=TIMENAME, FILENAME=os.path.join(DIRECTORY,local_file))
            # set permissions mode to MODE
            os.chmod(os.path.join(DIRECTORY,local_file), MODE)

    # close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

# PURPOSE: get attributes for a variable
def ncdf_attributes(fileID, var):
    # dictionary of attributes and list of attributes to retrieve
    attributes = {}
    attributes_list = ['calendar', 'long_name', 'positive',
        'standard_name', 'units', 'valid_range']
    # for each potential attribute
    for att_name in attributes_list:
        try:
            att_val = fileID[var].getncattr(att_name)
        except Exception as e:
            pass
        else:
            attributes[att_name] = att_val
    # return the dictionary of attributes
    return attributes

# PURPOSE: write output model layer fields data to file
def ncdf_model_write(dinput, attributes, fill_value,
    VARNAME=None, TNAME=None, QNAME=None, LONNAME=None, LATNAME=None,
    LEVELNAME=None, TIMENAME=None, FILENAME=None):
    # opening NetCDF4 file for writing
    fileID = netCDF4.Dataset(FILENAME, 'w', format="NETCDF4")

    # Defining the NetCDF4 dimensions and creating dimension variables
    nc = {}
    for key in [LONNAME,LATNAME,TIMENAME,LEVELNAME]:
        fileID.createDimension(key, len(dinput[key]))
        nc[key] = fileID.createVariable(key,dinput[key].dtype,(key,))
    # creating the layered NetCDF4 variables
    for key in [TNAME,QNAME]:
        nc[key] = fileID.createVariable(key, dinput[key].dtype,
            (TIMENAME,LEVELNAME,LATNAME,LONNAME,),
            fill_value=fill_value, zlib=True)
    # creating the surface NetCDF4 variables
    for key in [VARNAME]:
        nc[key] = fileID.createVariable(key, dinput[key].dtype,
            (TIMENAME,LATNAME,LONNAME,),
            fill_value=fill_value, zlib=True)

    # filling NetCDF4 variables
    for key,val in dinput.items():
        nc[key][:] = val.copy()
        # set netCDF4 attributes for variable
        for att_name,att_val in attributes[key].items():
            nc[key].setncattr(att_name, att_val)

    # Defining file-level attributes
    for att_name,att_val in attributes['ROOT'].items():
        fileID.setncattr(att_name, att_val)
    # add software information
    fileID.software_reference = mdlhmc.version.project_name
    fileID.software_version = mdlhmc.version.full_version
    # date created
    fileID.date_created = time.strftime('%Y-%m-%d',time.localtime())

    # Output NetCDF structure information
    logging.info(os.path.basename(FILENAME))
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
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # MERRA-2 product shortname
    parser.add_argument('--shortname','-s',
        type=str, default='M2I3NVASM',
        help='MERRA-2 product shortname')
    # MERRA-2 version
    parser.add_argument('--version','-v',
        type=str, default='5.12.4',
        help='MERRA-2 version')
    # years to download
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2000,now.tm_year+1),
        help='Years of model outputs to sync')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # NASA_GESDISC_MERRA2_monthly_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    # get NASA Earthdata credentials
    try:
        args.user,_,args.password = netrc.netrc(args.netrc).authenticators(URS)
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
    mdlhmc.utilities.build_opener(args.user, args.password,
        password_manager=True, authorization_header=False)

    # check internet connection before attempting to run program
    HOST = 'https://goldsmr5.gesdisc.eosdis.nasa.gov/'
    if mdlhmc.utilities.check_credentials(HOST):
        gesdisc_merra_monthly(args.directory, args.shortname,
            VERSION=args.version,
            YEARS=args.year,
            TIMEOUT=args.timeout,
            LOG=args.log,
            VERBOSE=args.verbose,
            MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
