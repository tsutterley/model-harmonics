#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (10/2024)
Download and management utilities for syncing time and auxiliary files
Adds additional modules to the gravity_toolkit utilities

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python (https://pypi.python.org/pypi/lxml)
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2024: update CMR search utility to replace deprecated scrolling
        https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html
    Updated 11/2023: updated ssl context to fix deprecation error
    Updated 05/2023: use pathlib to define and operate on paths
    Updated 01/2023: add default ssl context attribute with protocol
    Updated 12/2022: functions for managing and maintaining git repositories
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added option to use CMR provided GES DISC subsetting host
    Updated 08/2022: hardcode GES DISC subsetting API hostname
    Updated 06/2022: add NASA Common Metadata Repository (CMR) queries
        added function to build GES DISC subsetting API requests
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 01/2021
"""
# extend gravity_toolkit utilities
from __future__ import annotations
from gravity_toolkit.utilities import *

# PURPOSE: get the git hash value
def get_git_revision_hash(
        refname: str = 'HEAD',
        short: bool = False
    ):
    """
    Get the ``git`` hash value for a particular reference

    Parameters
    ----------
    refname: str, default HEAD
        Symbolic reference name
    short: bool, default False
        Return the shorted hash value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'rev-parse']
    cmd.append('--short') if short else None
    cmd.append(refname)
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

# PURPOSE: get the current git status
def get_git_status():
    """Get the status of a ``git`` repository as a boolean value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'status', '--porcelain']
    with warnings.catch_warnings():
        return bool(subprocess.check_output(cmd))

def _create_default_ssl_context() -> ssl.SSLContext:
    """Creates the default SSL context
    """
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    _set_ssl_context_options(context)
    context.options |= ssl.OP_NO_COMPRESSION
    return context

def _create_ssl_context_no_verify() -> ssl.SSLContext:
    """Creates an SSL context for unverified connections
    """
    context = _create_default_ssl_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    return context

def _set_ssl_context_options(context: ssl.SSLContext) -> None:
    """Sets the default options for the SSL context
    """
    if sys.version_info >= (3, 10) or ssl.OPENSSL_VERSION_INFO >= (1, 1, 0, 7):
        context.minimum_version = ssl.TLSVersion.TLSv1_2
    else:
        context.options |= ssl.OP_NO_SSLv2
        context.options |= ssl.OP_NO_SSLv3
        context.options |= ssl.OP_NO_TLSv1
        context.options |= ssl.OP_NO_TLSv1_1

# default ssl context
_default_ssl_context = _create_ssl_context_no_verify()

# PURPOSE: list a directory on NASA GES DISC https server
def gesdisc_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        build: bool = False,
        timeout: int | None = None,
        urs: str = 'urs.earthdata.nasa.gov',
        parser = lxml.etree.HTMLParser(),
        format: str = '%Y-%m-%d %H:%M',
        pattern: str = '',
        sort: bool = False
    ):
    """
    List a directory on NASA GES DISC servers

    Parameters
    ----------
    HOST: str or list
        remote https host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener with NASA Earthdata credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        Earthdata login URS 3 host
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
    format: str, default '%Y-%m-%d %H:%M'
        format for input time string
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    # build urllib2 opener with credentials
    if build:
        build_opener(username, password, password_manager=True,
            authorization_header=False)
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from https
    try:
        # Create and submit request.
        request=urllib2.Request(posixpath.join(*HOST))
        response=urllib2.urlopen(request, timeout=timeout)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response,parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        # get the Unix timestamp value for a modification time
        lastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            lastmod = [lastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames,lastmod)

# PURPOSE: filter the CMR json response for desired data files
def cmr_filter_json(
        search_results: dict,
        endpoint: str = "data",
        request_type: str = "application/x-netcdf"
    ):
    """
    Filter the CMR json response for desired data files

    Parameters
    ----------
    search_results: dict
        json response from CMR query
    endpoint: str, default 'data'
        url endpoint type

            - ``'data'``: NASA Earthdata https archive
            - ``'opendap'``: NASA Earthdata OPeNDAP archive
            - ``'s3'``: NASA Earthdata Cumulus AWS S3 bucket
    request_type: str, default 'application/x-netcdf'
        data type for reducing CMR query

    Returns
    -------
    granule_names: list
        Model granule names
    granule_urls: list
        Model granule urls
    granule_mtimes: list
        Model granule modification times
    """
    # output list of granule ids, urls and modified times
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    # check that there are urls for request
    if ('feed' not in search_results) or ('entry' not in search_results['feed']):
        return (granule_names,granule_urls)
    # descriptor links for each endpoint
    rel = {}
    rel['data'] = "http://esipfed.org/ns/fedsearch/1.1/data#"
    rel['opendap'] = "http://esipfed.org/ns/fedsearch/1.1/service#"
    rel['s3'] = "http://esipfed.org/ns/fedsearch/1.1/s3#"
    # iterate over references and get cmr location
    for entry in search_results['feed']['entry']:
        granule_names.append(entry['producer_granule_id'])
        granule_mtimes.append(get_unix_time(entry['updated'],
            format='%Y-%m-%dT%H:%M:%S.%f%z'))
        for link in entry['links']:
            # skip inherited granules
            if ('inherited' in link.keys()):
                continue
            # append if selected endpoint
            if (link['rel'] == rel[endpoint]):
                granule_urls.append(link['href'])
                break
            # alternatively append if selected data type
            if ('type' not in link.keys()):
                continue
            if (link['type'] == request_type):
                granule_urls.append(link['href'])
                break
    # return the list of urls, granule ids and modified times
    return (granule_names, granule_urls, granule_mtimes)

# PURPOSE: cmr queries for model data products
def cmr(
        short_name: str,
        version: str | int | None = None,
        start_date: str | None = None,
        end_date: str | None = None,
        provider: str | None = 'GES_DISC',
        endpoint: str | None = 'data',
        request_type: str | None = 'application/x-netcdf',
        verbose: bool = False,
        fid = sys.stdout
    ):
    """
    Query the NASA Common Metadata Repository (CMR) for model data

    Parameters
    ----------
    short_name: str
        Model shortname in the CMR system
    version: str or NoneType, default None
        Model version
    start_date: str or NoneType, default None
        starting date for CMR product query
    end_date: str or NoneType, default None
        ending date for CMR product query
    provider: str, default 'GES_DISC'
        CMR data provider

            - ``'GES_DISC'``: GESDISC
            - ``'GESDISCCLD'``: GESDISC Cumulus
            - ``'PODAAC'``: PO.DAAC Drive
            - ``'POCLOUD'``: PO.DAAC Cumulus
    endpoint: str, default 'data'
        url endpoint type

            - ``'data'``: NASA Earthdata https archive
            - ``'opendap'``: NASA Earthdata OPeNDAP archive
            - ``'s3'``: NASA Earthdata Cumulus AWS S3 bucket
    request_type: str, default 'application/x-netcdf'
        data type for reducing CMR query
    verbose: bool, default False
        print CMR query information
    fid: obj, default sys.stdout
        open file object to print if verbose

    Returns
    -------
    granule_names: list
        Model granule names
    granule_urls: list
        Model granule urls
    granule_mtimes: list
        Model granule modification times
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # build urllib2 opener with SSL context
    # https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    # Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    handler.append(urllib2.HTTPSHandler(context=_default_ssl_context))
    # create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    # build CMR query
    cmr_query_type = 'granules'
    cmr_format = 'json'
    cmr_page_size = 2000
    CMR_HOST = ['https://cmr.earthdata.nasa.gov','search',
        f'{cmr_query_type}.{cmr_format}']
    # build list of CMR query parameters
    CMR_KEYS = []
    CMR_KEYS.append(f'?provider={provider}')
    CMR_KEYS.append('&sort_key[]=start_date')
    CMR_KEYS.append('&sort_key[]=producer_granule_id')
    CMR_KEYS.append(f'&page_size={cmr_page_size}')
    # dictionary of product shortnames and version
    CMR_KEYS.append(f'&short_name={short_name}')
    if version:
        CMR_KEYS.append(f'&version={version}')
    # append keys for start and end time
    # verify that start and end times are in ISO format
    start_date = isoformat(start_date) if start_date else ''
    end_date = isoformat(end_date) if end_date else ''
    CMR_KEYS.append(f'&temporal={start_date},{end_date}')
    # full CMR query url
    cmr_query_url = "".join([posixpath.join(*CMR_HOST),*CMR_KEYS])
    logging.info(f'CMR request={cmr_query_url}')
    # output list of granule names and urls
    granule_names = []
    granule_urls = []
    granule_mtimes = []
    cmr_search_after = None
    while True:
        req = urllib2.Request(cmr_query_url)
        # add CMR search after header
        if cmr_search_after:
            req.add_header('CMR-Search-After', cmr_search_after)
            logging.debug(f'CMR-Search-After: {cmr_search_after}')
        response = opener.open(req)
        # get search after index for next iteration
        headers = {k.lower():v for k,v in dict(response.info()).items()}
        cmr_search_after = headers.get('cmr-search-after')
        # read the CMR search as JSON
        search_page = json.loads(response.read().decode('utf8'))
        ids,urls,mtimes = cmr_filter_json(search_page,
            endpoint=endpoint, request_type=request_type)
        if not urls or cmr_search_after is None:
            break
        # extend lists
        granule_names.extend(ids)
        granule_urls.extend(urls)
        granule_mtimes.extend(mtimes)
    # return the list of granule ids, urls and modification times
    return (granule_names, granule_urls, granule_mtimes)

# PURPOSE: build requests for the GES DISC subsetting API
def build_request(
        short_name: str,
        dataset_version: str | int,
        url: str | None,
        host: str | None = None,
        variables: list | None = [],
        format: str | None = 'bmM0Lw',
        service: str | None = 'L34RS_MERRA2',
        version: str | None = '1.02',
        bbox: list | None = [-90,-180,90,180],
        **kwargs
    ):
    """
    Build requests for the GES DISC subsetting API

    Parameters
    ----------
    short_name: str
        Model shortname in the CMR system
    dataset_version: str
        Model version
    url: str
        url for granule returned by the CMR system
    host: str or NoneType, default None
        Override host provider for GES DISC subsetting

        Default is host provider given by CMR request
    variables: list, default []
        Variables for product to subset
    format: str, default 'bmM0Lw'
        Coded output format for GES DISC subsetting API
    service: str, default 'L34RS_MERRA2'
        GES DISC subsetting API service
    version: str, default '1.02'
        GES DISC subsetting API service version
    bbox: list, default [-90,-180,90,180]
        Bounding box to spatially subset
    **kwargs: dict, default {}
        Additional parameters for GES DISC subsetting API

    Returns
    -------
    request_url: str
        Formatted url for GES DISC subsetting API
    """
    # split CMR supplied url for granule
    HOST,*args = url_split(url)
    host = HOST if (host is None) else host
    api_host = posixpath.join(host,'daac-bin','OTF','HTTP_services.cgi?')
    # create parameters to be encoded
    kwargs['FILENAME'] = posixpath.join(posixpath.sep, *args)
    kwargs['FORMAT'] = format
    kwargs['SERVICE'] = service
    kwargs['VERSION'] = version
    kwargs['BBOX'] = ','.join(map(str, bbox))
    kwargs['SHORTNAME'] = short_name
    kwargs['DATASET_VERSION'] = dataset_version
    kwargs['VARIABLES'] = ','.join(variables)
    # return the formatted request url
    request_url = api_host + urlencode(kwargs)
    return request_url
