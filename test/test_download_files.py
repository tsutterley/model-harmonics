#!/usr/bin/env python
"""
test_download_and_read.py (08/2020)
Tests that files can be downloaded from data sources
"""

import os
import pytest
import inspect
import model_harmonics as mdlhmc

filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))


# PURPOSE: Download an ECCO realtime file from JPL ECCO Drive
# ECCO-KFS data might also be deprecated and unavailable for download
@pytest.mark.skip(reason='Need to update to PO.DAAC Cumulus')
def test_ECCO_realtime_download(username, webdav):
    HOST = [
        'https://ecco.jpl.nasa.gov',
        'drive',
        'files',
        'NearRealTime',
        'KalmanFilter',
        'kf080i_1993',
        'n10day_01_09',
        'OBP_08_08.00001_02160_012.cdf',
    ]
    local = os.path.join(filepath, HOST[-1])
    try:
        # build opener for JPL ECCO Drive
        mdlhmc.utilities.build_opener(username, webdav)
        # download ECCO file
        mdlhmc.utilities.from_drive(
            HOST,
            build=False,
            verbose=True,
            local=local,
        )
    except mdlhmc.utilities.urllib2.URLError as exc:
        pytest.xfail(exc.reason)
    except mdlhmc.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    except ValueError as exc:
        pytest.xfail(exc)
    # check that the file was downloaded
    assert os.access(local, os.F_OK)
    # clean up files
    os.remove(local)


# PURPOSE: Download a GLDAS file from NASA GESDISC
def test_GLDAS_NOAH10M_download(username, password):
    HOST = [
        'https://hydro1.gesdisc.eosdis.nasa.gov',
        'data',
        'GLDAS',
        'GLDAS_NOAH10_M.2.1',
        '2000',
        'GLDAS_NOAH10_M.A200001.021.nc4',
    ]
    local = os.path.join(filepath, HOST[-1])
    # try to download GLDAS file from NASA GESDISC
    try:
        # build opener for GESDISC
        mdlhmc.utilities.build_opener(
            username,
            password,
            password_manager=True,
            authorization_header=False,
        )
        # download GLDAS file
        mdlhmc.utilities.from_http(
            HOST,
            context=None,
            verbose=True,
            local=local,
        )
    except mdlhmc.utilities.urllib2.URLError as exc:
        pytest.xfail(exc.reason)
    except mdlhmc.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    except ValueError as exc:
        pytest.xfail(exc)
    # check that the file was downloaded
    assert os.access(local, os.F_OK)
    # clean up files
    os.remove(local)


# PURPOSE: Download a MERRA-2 SMB file from NASA GESDISC
def test_MERRA2_SMB_download(username, password):
    HOST = [
        'https://goldsmr4.gesdisc.eosdis.nasa.gov',
        'data',
        'MERRA2_MONTHLY',
        'M2TMNXGLC.5.12.4',
        '1980',
        'MERRA2_100.tavgM_2d_glc_Nx.198001.nc4',
    ]
    local = os.path.join(filepath, HOST[-1])
    # try to download MERRA-2 file from NASA GESDISC
    try:
        # build opener for GESDISC
        mdlhmc.utilities.build_opener(
            username,
            password,
            password_manager=True,
            authorization_header=False,
        )
        # download MERRA-2 file
        mdlhmc.utilities.from_http(
            HOST,
            context=None,
            verbose=True,
            local=local,
        )
    except mdlhmc.utilities.urllib2.URLError as exc:
        pytest.xfail(exc.reason)
    except mdlhmc.utilities.urllib2.HTTPError as exc:
        pytest.xfail(exc.reason)
    except ValueError as exc:
        pytest.xfail(exc)
    # check that the file was downloaded
    assert os.access(local, os.F_OK)
    # clean up files
    os.remove(local)
