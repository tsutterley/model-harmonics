#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (01/2021)
Download and management utilities for syncing time and auxiliary files
Adds additional modules to the gravity_toolkit utilities

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python (https://pypi.python.org/pypi/lxml)
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Written 01/2021
"""
#-- extend gravity_toolkit utilities
from gravity_toolkit.utilities import *

#-- PURPOSE: list a directory on NASA GES DISC https server
def gesdisc_list(HOST,username=None,password=None,build=False,timeout=None,
    urs='urs.earthdata.nasa.gov',parser=lxml.etree.HTMLParser(),
    format='%Y-%m-%d %H:%M',pattern='',sort=False):
    """
    List a directory on NASA GES DISC servers

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    username: NASA Earthdata username
    password: NASA Earthdata password
    build: Build opener with NASA Earthdata credentials
    timeout: timeout in seconds for blocking operations
    urs: Earthdata login URS 3 host
    parser: HTML parser for lxml
    format: format for input time string
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    colnames: list of column names in a directory
    collastmod: list of last modification times for items in the directory
    """
    #-- use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener with credentials
    if build:
        build_opener(username, password, password_manager=True,
            authorization_header=False)
    #-- try listing from https
    try:
        #-- Create and submit request.
        request=urllib2.Request(posixpath.join(*HOST))
        response=urllib2.urlopen(request,timeout=timeout)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response,parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        #-- get the Unix timestamp value for a modification time
        lastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            #-- reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            #-- sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        #-- return the list of column names and last modified times
        return (colnames,lastmod)
