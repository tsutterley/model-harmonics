#!/usr/bin/env python
u"""
jpl_ecco_webdav.py
Written by Tyler Sutterley (12/2022)

Retrieves and prints a user's JPL ECCO Drive WebDAV credentials

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add ECCO Drive to NASA Earthdata Applications and get WebDAV Password
https://ecco.jpl.nasa.gov/drive

CALLING SEQUENCE:
    python jpl_ecco_webdav.py --user <username>
    where <username> is your NASA Earthdata credentials

OUTPUTS:
    JPL ECCO Drive WebDAV credentials

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for authentication
    -A, --append: append .netrc file instead of printing

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
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 05/2021: use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Written 12/2020 for public release
"""
from __future__ import print_function

import sys
import os
import netrc
import base64
import getpass
import builtins
import argparse
import posixpath
import lxml.etree
import gravity_toolkit as gravtk

# PURPOSE: retrieve JPL ECCO Drive WebDAV credentials
def jpl_ecco_webdav(USER, PASSWORD, parser=lxml.etree.HTMLParser()):
    # build opener for retrieving JPL ECCO Drive WebDAV credentials
    # Add the username and password for NASA Earthdata Login system
    URS = 'https://urs.earthdata.nasa.gov'
    gravtk.utilities.build_opener(USER, PASSWORD,
        password_manager=True, authorization_header=True, urs=URS)
    # All calls to urllib2.urlopen will now use handler
    # Make sure not to include the protocol in with the URL, or
    # HTTPPasswordMgrWithDefaultRealm will be confused.
    HOST = posixpath.join('https://ecco.jpl.nasa.gov','drive')
    parameters = gravtk.utilities.urlencode(
        {'client_id':'gA5gkD03X2RMcpJUi8zbRA', 'response_type':'code',
        'state':base64.b64encode(HOST.encode()),
        'redirect_uri':posixpath.join(HOST,'authenticated'),
        'required_scope': 'country+study_area'}
    )
    # retrieve cookies from NASA Earthdata URS
    request = gravtk.utilities.urllib2.Request(
        url=posixpath.join(URS,'oauth',f'authorize?{parameters}'))
    gravtk.utilities.urllib2.urlopen(request)
    # read and parse request for webdav password
    request = gravtk.utilities.urllib2.Request(url=HOST)
    response = gravtk.utilities.urllib2.urlopen(request,timeout=20)
    tree = lxml.etree.parse(response, parser)
    WEBDAV, = tree.xpath('//input[@id="password"]/@value')
    # return webdav password
    return WEBDAV

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Retrieves and prints a user's JPL ECCO WebDAV
            credentials
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
    # append to netrc
    parser.add_argument('--append','-A',
        default=False, action='store_true',
        help='Append .netrc file instead of printing')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # NASA Earthdata hostname
    URS = 'urs.earthdata.nasa.gov'
    # JPL JPL ECCO drive hostname
    HOST = 'ecco.jpl.nasa.gov'
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

    # check internet connection before attempting to run program
    DRIVE = posixpath.join('https://ecco.jpl.nasa.gov','drive')
    if gravtk.utilities.check_connection(DRIVE):
        # compile HTML parser for lxml
        WEBDAV = jpl_ecco_webdav(args.user, args.password)
        # output to terminal or append to netrc file
        if args.append:
            # append to netrc file and set permissions level
            with open(args.netrc, mode='a+') as f:
                f.write(f'machine {args.user} login {HOST} password {WEBDAV}\n')
                os.chmod(args.netrc, 0o600)
        else:
            print(f'\nWebDAV Password for {args.user}@{HOST}:\n\t{WEBDAV}')

# run main program
if __name__ == '__main__':
    main()
