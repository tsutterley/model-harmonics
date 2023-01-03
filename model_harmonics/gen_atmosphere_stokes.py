#!/usr/bin/env python
u"""
gen_atmosphere_stokes.py
Written by Tyler Sutterley (01/2023)
Calculates spherical harmonic fields from 3D atmospheric geopotential
    height and pressure difference fields

CALLING SEQUENCE:
    Ylms = gen_atmosphere_stokes(GPH, pressure, lon, lat, LMAX=60,
        ELLIPSOID=ELLIPSOID, GEOID=GEOID, PLM=PLM, LOVE=(hl,kl,ll))

INPUTS:
    GPH: geopotential heights at model levels
    pressure: pressure differences between model levels
    lon: longitude array
    lat: latitude array

OUTPUTS:
    Ylms: harmonics object
        clm: fully-normalized cosine spherical harmonic coefficients
        slm: fully-normalied sine spherical harmonic coefficients
        l: spherical harmonic degree to LMAX
        m: spherical harmonic order to MMAX

OPTIONS:
    LMAX: Upper bound of Spherical Harmonic Degrees (default = 60)
    MMAX: Upper bound of Spherical Harmonic Orders (default = LMAX)
    ELLIPSOID: reference ellipsoid name
    GEOID: geoid height
    PLM: input Legendre polynomials
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)
    METHOD: method of integrating over pressure levels
        SW02: Swenson and Wahr (2002)
        BC05: Boy and Chao (2005) (default)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units
    constants.py: calculate reference parameters for common ellipsoids
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCE:
    JP Boy and B Chao, Precise evaluation of atmospheric loading effects on
    Earth's time-variable gravity field, Journal of Geophysical Research:
    Solid Earth, 110(B8), 2005. https://doi.org/10.1029/2002JB002333

    S Swenson and J Wahr, Estimated effects of the vertical structure of
    atmospheric mass on the time-variable geoid, Journal of Geophysical
    Research: Solid Earth, 107(B9), 2002. https://doi.org/10.1029/2000JB000024

    S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order
    normalised associated Legendre functions" Journal of Geodesy,
    76: 279-299, 2002. https://doi.org/10.1007/s00190-002-0216-2

UPDATE HISTORY:
    Updated 01/2023: refactored associated legendre polynomials
    Updated 12/2022: constants class in place of geoid-toolkit ref_ellipsoid
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 01/2021: added function docstrings
    Updated 05/2020: use harmonics class for spherical harmonic operations
    Updated 04/2020: made Legendre polynomials and Love numbers options
        using the units class for converting to normalized spherical harmonics
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Written 03/2018
"""

import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes
from model_harmonics.constants import constants

# PURPOSE: calculates spherical harmonic fields from 3D atmospheric
# geopotential height and pressure difference fields
def gen_atmosphere_stokes(GPH, pressure, lon, lat, LMAX=60, MMAX=None,
    ELLIPSOID=None, GEOID=None, PLM=None, LOVE=None, METHOD='BC05'):
    """
    Converts 3D atmospheric geopotential height and pressure difference
    fields from the spatial domain to spherical harmonic coefficients

    Parameters
    ----------
    GPH: float
        geopotential heights at model levels
    pressure: float
        pressure differences between model levels
    lon: float
        longitude array
    lat: float
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    ELLIPSOID: str or NoneType, default None
        reference ellipsoid name
    GEOID: float or NoneType, default None
        geoid height
    PLM: float or NoneType, default None
        Legendre polynomials
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    METHOD: str, default 'BC05'
        Method of integrating over pressure levels

            - ``'BC05'``: [Boy2005]_
            - ``'SW02'``: [Swenson2002]_

    Returns
    -------
    clm: float
        fully-normalized cosine spherical harmonic coefficients
    slm: float
        fully-normalized sine spherical harmonic coefficients
    l: int
        spherical harmonic degree to LMAX
    m: int
        spherical harmonic order to MMAX

    References
    ----------
    .. [Boy2005] J.-P. Boy and B. F. Chao, "Precise evaluation of
        atmospheric loading effects on Earth's time‐variable gravity field",
        *Journal of Geophysical Research: Solid Earth*, 110(B08412), (2005).
        `doi: 10.1029/2002JB002333 <https://doi.org/10.1029/2002JB002333>`_
    .. [Swenson2002] S. Swenson and J. Wahr, "Estimated effects of the vertical
        structure of atmospheric mass on the time‐variable geoid",
        *Journal of Geophysical Research*, 107(B9), 2194, (2002).
        `doi: 10.1029/2000JB000024 <https://doi.org/10.1029/2000JB000024>`_
    """

    # converting LMAX to integer
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX

    # number of pressure levels, longitudes and latitudes
    nlevels,nlat,nlon = np.shape(GPH)
    # grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    # longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    # colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    # calculate longitudes and colatitudes in radians
    phi = lon*np.pi/180.0
    phi = np.squeeze(phi)[np.newaxis,:]
    th = (90.0 - np.squeeze(lat))*np.pi/180.0
    # calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    # Earth Parameters
    ellipsoid_params = constants(ellipsoid=ELLIPSOID)
    # semimajor axis of ellipsoid [m]
    a_axis = ellipsoid_params.a_axis
    # ellipsoidal flattening
    flat = ellipsoid_params.flat
    # Average Radius of the Earth having the same volume [m]
    rad_e = ellipsoid_params.rad_e
    # first numerical eccentricity
    ecc1 = ellipsoid_params.ecc1
    # convert from geodetic latitude to geocentric latitude
    # prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.cos(gridtheta)**2.0)

    # Coefficient for calculating Stokes coefficients from pressure field
    # SH Degree dependent factors with indirect loading components
    factors = gravity_toolkit.units(lmax=LMAX,a_axis=100.0*a_axis,flat=flat)
    dfactor = factors.spatial(*LOVE).mmwe
    # Multiplying sin(th) with differentials of theta and phi
    # to calculate the integration factor at each latitude
    int_fact = np.sin(th)*dphi*dth

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    m = np.arange(MMAX+1)
    ccos = np.cos(np.dot(m[:, np.newaxis],phi))
    ssin = np.sin(np.dot(m[:, np.newaxis],phi))

    # added option to precompute plms to improve computational speed
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # Fully-normalized Legendre Polynomials
    # Multiplying by the units conversion factor (conv) to
    # Multiplying by integration factors [sin(theta)*dtheta*dphi]
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    for j in range(0,nlat):
        plm[:,m,j] = PLM[:,m,j]*int_fact[j]

    # gravitational acceleration at the Earth's mean spherical surface
    g0 = 9.80665
    # gravitational acceleration at the equator and at mean sea level
    ge = 9.780356
    # gravitational acceleration at the mean sea level over gridtheta
    gs = ge*(1.0+5.2885e-3*np.cos(gridtheta)**2-5.9e-6*np.cos(2.0*gridtheta)**2)
    # calculate radii and gravity for each pressure level
    R = np.zeros((nlevels,nlat,nlon))
    gamma_h = np.zeros((nlevels,nlat,nlon))
    for p in range(nlevels):
        # orthometric height from List (1958)
        # as described in Boy and Chao (2005)
        orthometric = (1.0 - 0.002644*np.cos(2.0*gridtheta))*GPH[p,:,:] + \
            (1.0 - 0.0089*np.cos(2.0*gridtheta))*(GPH[p,:,:]**2)/6.245e6
        # calculate X, Y and Z from geodetic latitude and longitude
        X = (N + GEOID + orthometric) * np.sin(gridtheta) * np.cos(gridphi)
        Y = (N + GEOID + orthometric) * np.sin(gridtheta) * np.sin(gridphi)
        Z = (N * (1.0 - ecc1**2.0) + GEOID + orthometric) * np.cos(gridtheta)
        # calculate radius of level
        R[p,:,:] = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
        # calculate normal gravity at each height above mean sea level
        gamma_h[p,:,:] = gs*(1.0-2.0*(1.006803-0.06706*np.cos(gridtheta)**2)*
            (orthometric/rad_e) + 3.0*(orthometric/rad_e)**2)

    # Initializing preliminary spherical harmonic matrices
    yclm = np.zeros((LMAX+1,MMAX+1))
    yslm = np.zeros((LMAX+1,MMAX+1))
    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    for l in range(0,LMAX+1):# equivalent to 0:LMAX
        mm = np.min([MMAX,l])# truncate to MMAX if specified (if l > MMAX)
        m = np.arange(0,mm+1)# mm+1 elements between 0 and mm
        # total pressure factor
        pfactor = np.zeros((nlat,nlon))
        # iterate over pressure levels
        for p in range(nlevels):
            # if using Swenson and Wahr (2002) or Boy and Chao (2005)
            if (METHOD == 'SW02'):
                # calculate pressure change/gravity ratio
                PG = pressure[p,:,:]/g0
                # add to pressure factor (pfactor) to integrate over levels
                pfactor += PG*(rad_e/(rad_e-GPH[p,:,:])+(GEOID/rad_e))**(l+4)
            elif (METHOD == 'BC05'):
                # calculate pressure change/gravity ratio
                PG = pressure[p,:,:]/gamma_h[p,:,:]
                # add to pressure factor (pfactor) to integrate over levels
                pfactor += PG*(R[p,:,:]/rad_e)**(l+2)
        # Multiplying gridded data with sin/cos of m#phis
        # This will sum through all phis in the dot product
        # need to reform pfactor to lonXlat as is originally latXlon
        # output [m,theta]
        dcos = np.dot(ccos,-np.transpose(pfactor))
        dsin = np.dot(ssin,-np.transpose(pfactor))
        # Summing product of plms and data over all latitudes
        # axis=1 signifies the direction of the summation (colatitude (th))
        # ycos and ysin are the SH coefficients before normalizing
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        # Multiplying by coefficients to normalize
        Ylms.clm[l,m] = dfactor[l]*yclm[l,m]
        Ylms.slm[l,m] = dfactor[l]*yslm[l,m]

    # return the harmonics object
    return Ylms
