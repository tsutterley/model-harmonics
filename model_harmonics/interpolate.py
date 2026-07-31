"""
interpolate.py
Written by Tyler Sutterley (07/2026)
Routines to interpolate data from pre-computed spatial grids

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/

UPDATE HISTORY:
    Written 07/2026
"""

from __future__ import annotations

import logging
import numpy as np
import netCDF4
from geoid_toolkit.interpolate import Interpolate
from model_harmonics.utilities import import_dependency, get_data_path

# attempt imports
netCDF4 = import_dependency('netCDF4')


# PURPOSE: read ocean depth netCDF4 file and interpolate to points
def ocean_depth(
    longitude: np.ndarray,
    latitude: np.ndarray,
    model: str = '2020',
    resolution: str = '1440x720',
    method: str = 'linear',
    **kwargs,
):
    """
    Interpolate ocean depth to input lat/lon points

    Parameters
    ----------
    longitude: np.ndarray
        Longitudes of output points
    latitude: np.ndarray
        Latitudes of output points
    model: str, default "2020"
        GEBCO bathymetry model to use in the interpolation
    resolution: str, default "1440x720"
        Spatial resolution of the bathymetry model
    method: str, default "linear"
        Interpolation method to use

        - ``'linear'``: linear interpolation
        - ``'nearest'``: nearest-neighbor interpolation
    kwargs: keyword arguments
        Additional keyword arguments for the ``Interpolate`` class

    Returns
    -------
    depth: np.ndarray
        Ocean depth at input lat/lon points
    """
    # read bathymetry data
    FILENAME = get_data_path(['data', f'DEPTH.{model}.{resolution}.nc'])
    logging.debug(str(FILENAME))
    with netCDF4.Dataset(FILENAME, mode='r') as fileID:
        z = fileID.variables['depth'][:].copy()
        lon = fileID.variables['lon'][:].copy()
        lat = fileID.variables['lat'][:].copy()
    # interpolate values for a given method
    interp = Interpolate((lat, lon), z, method=method, **kwargs)
    # interpolate ocean depth for input lat/lon points
    depth = interp((latitude, longitude))
    # return ocean depth
    return depth
