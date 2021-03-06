#!/usr/bin/env python
u"""
gesdisc_merra_download.py
Written by Tyler Sutterley (04/2021)

This program downloads MERRA-2 products using a links list provided by the
    Goddard Earth Sciences Data and Information Server Center
    https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

Register with NASA Earthdata Login system:
    https://urs.earthdata.nasa.gov

Add "NASA GESDISC DATA ARCHIVE" to Earthdata Applications:
    https://urs.earthdata.nasa.gov/approve_app?client_id=e2WVk8Pw6weeLUKZYOxvTQ

CALLING SEQUENCE:
    python gesdisc_merra_download.py --user <username> links_list_file
    where <username> is your NASA Earthdata username

INPUTS:
    links_list_file: GES DISC generated file listing files to download

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -P X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: Working data directory
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
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 01/2021: use argparse to set command line parameters
        using utilities program to build opener
    Updated 09/2019: added ssl context to urlopen headers
    Updated 08/2019: new GESDISC server and links list file format
        increased timeout to 20 seconds
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
import argparse
import builtins
import gravity_toolkit.utilities

#-- PURPOSE: download MERRA-2 files from GESDISC using a links list file
def gesdisc_merra_download(base_dir, links_list_file, LOG=False,
    VERBOSE=False, MODE=None):
    #-- full path to MERRA-2 directory
    DIRECTORY = os.path.join(base_dir,'MERRA-2')
    #-- check if DIRECTORY exists and recursively create if not
    if (not os.access(os.path.join(DIRECTORY), os.F_OK)):
        os.makedirs(os.path.join(DIRECTORY), MODE)

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: NASA_GESDISC_MERRA2_download_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NASA_GESDISC_MERRA2_download_{0}.log'.format(today)
        fid = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('NASA MERRA-2 Sync Log ({0})'.format(today), file=fid)
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- read the links list file
    with open(links_list_file,'rb') as fileID:
        lines = fileID.read().decode("utf-8-sig").encode("utf-8").splitlines()

    #-- for each line in the links_list_file
    for f in lines:
        #-- extract filename from url
        HOST = gravity_toolkit.utilities.url_split(f.decode('utf-8'))
        if re.search(rb'LABEL\=(.*?)\&SHORTNAME',f):
            FILE, = re.findall(r'LABEL\=(.*?)\&SHORTNAME', f.decode('utf-8'))
        elif re.search(rb'MERRA2_(\d+)\.(.*?)\.(\d+)\.(.*?).nc',f):
            rx = re.compile(r'MERRA2_(\d+)\.(.*?)\.(\d+)\.(.*?).nc')
            MOD,DSET,YMD,AUX = rx.findall(f.decode('utf-8')).pop()
            FILE = 'MERRA2_{0}.{1}.{2}.SUB.nc'.format(MOD,DSET,YMD)
        else:
            FILE = HOST[-1]
        #-- output local file
        local_file = os.path.join(DIRECTORY,FILE)
        MD5 = gravity_toolkit.utilities.get_hash(local_file)
        #-- Create and submit request. There are a wide range of exceptions
        #-- that can be thrown here, including HTTPError and URLError.
        gravity_toolkit.utilities.from_http(HOST, context=None,
            local=local_file, hash=MD5, fid=fid, verbose=VERBOSE,
            mode=MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- Main program that calls gesdisc_merra_download()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Downloads MERRA-2 products using a links list
            provided by the Goddard Earth Sciences Data and Information
            Server Center (GES DISC)
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='GESDISC links list file')
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-P',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Output log file in form
    #-- NASA_GESDISC_MERRA2_download_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- print information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    #-- get NASA Earthdata credentials
    if not args.user and not os.access(args.netrc,os.F_OK):
        #-- check that NASA Earthdata credentials were entered
        args.user = builtins.input('Username for {0}: '.format(URS))
        #-- enter password securely from command-line
        args.password=getpass.getpass('Password for {0}@{1}: '.format(args.user,URS))
    elif not args.user and os.access(args.netrc,os.F_OK):
        args.user,_,args.password=netrc.netrc(args.netrc).authenticators(URS)
    elif not args.password:
        #-- enter password securely from command-line
        args.password=getpass.getpass('Password for {0}@{1}: '.format(args.user,URS))

    #-- build a urllib opener for NASA GESDISC
    #-- Add the username and password for NASA Earthdata Login system
    gravity_toolkit.utilities.build_opener(args.user, args.password,
        password_manager=True, authorization_header=False)

    #-- check internet connection before attempting to run program
    HOST = 'http://disc.sci.gsfc.nasa.gov/'
    if gravity_toolkit.utilities.check_connection(HOST):
        #-- for each links list file from GESDISC
        for FILE in args.file:
            gesdisc_merra_download(args.directory, FILE, LOG=args.log,
                VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
