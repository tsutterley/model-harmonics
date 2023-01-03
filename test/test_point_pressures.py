#!/usr/bin/env python
u"""
test_point_pressures.py (04/2021)
"""
import pytest
import numpy as np
import gravity_toolkit.units
from gravity_toolkit.utilities import get_data_path
from gravity_toolkit.read_love_numbers import read_love_numbers
from model_harmonics.gen_point_pressure import gen_point_pressure
from model_harmonics.gen_pressure_stokes import gen_pressure_stokes

# parameterize the number of point masses
@pytest.mark.parametrize("NPTS", np.random.randint(2,2000,size=1))
def test_point_masses(NPTS):
    # create spatial grid
    dlon,dlat = (1.0,1.0)
    lat = np.arange(90.0 - dlat/2.0, -90.0 - dlat/2.0, -dlat)
    lon = np.arange(-180.0 + dlon/2.0, 180.0 + dlon/2.0, dlon)
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridtheta = (90.0 - gridlat)*np.pi/180.0
    nlat,nlon = np.shape(gridlon)
    # longitude and colatitude degree spacings in radians
    dphi = dlon*np.pi/180.0
    dth = dlat*np.pi/180.0
    # radius at each point in meters (using a spherical geometry)
    rad_e = np.atleast_1d(gravity_toolkit.units().rad_e/100.0)

    # parameterize point masses
    LAT = lat[0]-dlat*np.random.randint(0,nlat,size=NPTS)
    LON  = lon[0]+dlon*np.random.randint(0,nlon,size=NPTS)
    th = (90.0 - LAT)*np.pi/180.0
    PRESSURE = 100.0 - 200.0*np.random.randn(NPTS)
    # calculate areas
    AREA = (rad_e**2)*np.sin(th)*dphi*dth

    # create test gridded field
    data = np.zeros((nlat,nlon))
    for i in range(NPTS):
        indy,indx = np.nonzero((gridlat == LAT[i]) & (gridlon == LON[i]))
        data[indy,indx] += PRESSURE[i]
    # gravitational acceleration at the equator and at mean sea level
    ge = 9.780356
    # gravitational acceleration at the mean sea level over thetas
    G = ge*(1.0+5.2885e-3*np.cos(th)**2-5.9e-6*np.cos(2.0*th)**2)
    gs = ge*(1.0+5.2885e-3*np.cos(gridtheta)**2-5.9e-6*np.cos(2.0*gridtheta)**2)

    # path to load Love numbers file
    love_numbers_file = get_data_path(['data','love_numbers'])
    # read load Love numbers
    hl,kl,ll = read_love_numbers(love_numbers_file)
    # calculate harmonics and degree amplitudes for each case
    grid_Ylms = gen_pressure_stokes(data, gs, rad_e, lon, lat,
        LMAX=60, LOVE=(hl,kl,ll))
    point_Ylms = gen_point_pressure(PRESSURE, G, rad_e, LON, LAT,
        AREA=AREA, LMAX=60, LOVE=(hl,kl,ll))

    # check that harmonic data is equal to machine precision
    difference_Ylms = grid_Ylms.copy()
    difference_Ylms.subtract(point_Ylms)
    harmonic_eps = np.finfo(np.float32).eps
    assert np.all(np.abs(difference_Ylms.clm) < harmonic_eps)
    # verify that the degree amplitudes are within tolerance
    assert np.all(np.abs(grid_Ylms.amplitude - point_Ylms.amplitude) < harmonic_eps)
