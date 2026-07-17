#!/usr/bin/env python
"""
gen_atmosphere_stokes.py
Written by Tyler Sutterley (07/2026)
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
    WEIGHT: custom latitudinal weighting function for gridded data
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
    datum.py: calculate reference parameters for common ellipsoids
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
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
    Updated 03/2023: improve typing for variables in docstrings
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
import gravity_toolkit as gravtk
from model_harmonics.datum import datum


# PURPOSE: calculates spherical harmonic fields from 3D atmospheric
# geopotential height and pressure difference fields
def gen_atmosphere_stokes(
    GPH,
    pressure,
    lon,
    lat,
    LMAX=60,
    MMAX=None,
    ELLIPSOID=None,
    GEOID=None,
    WEIGHT=None,
    PLM=None,
    LOVE=None,
    METHOD='BC05',
):
    """
    Converts 3D atmospheric geopotential height and pressure difference
    fields from the spatial domain to spherical harmonic coefficients

    Parameters
    ----------
    GPH: np.ndarray
        geopotential heights at model levels
    pressure: np.ndarray
        pressure differences between model levels
    lon: np.ndarray
        longitude array
    lat: np.ndarray
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    ELLIPSOID: str or NoneType, default None
        Reference ellipsoid name
    GEOID: np.ndarray or NoneType, default None
        Geoid height
    WEIGHT: np.ndarray or NoneType, default None
        Custom latitudinal weighting function for gridded data
    PLM: np.ndarray or NoneType, default None
        Legendre polynomials
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)
    METHOD: str, default 'BC05'
        Method of integrating over pressure levels

            - ``'BC05'``: :cite:p:`Boy:2005el`
            - ``'SW02'``: :cite:p:`Swenson:2002kf`

    Returns
    -------
    clm: np.ndarray
        fully-normalized cosine spherical harmonic coefficients
    slm: np.ndarray
        fully-normalized sine spherical harmonic coefficients
    l: np.ndarray
        spherical harmonic degree to LMAX
    m: np.ndarray
        spherical harmonic order to MMAX
    """

    # converting LMAX to integer
    LMAX = np.int64(LMAX)
    # upper bound of spherical harmonic orders (default = LMAX)
    MMAX = np.copy(LMAX) if not MMAX else MMAX

    # number of pressure levels, longitudes and latitudes
    nlevels, nlat, nlon = np.shape(GPH)
    # calculate longitudes and colatitudes in radians
    phi = np.radians(np.squeeze(lon))
    th = np.radians(90.0 - np.squeeze(lat))
    # calculate meshgrid from latitude and colongitude
    gridphi, gridth = np.meshgrid(phi, th)

    # Earth Parameters
    ellipsoid_params = datum(ellipsoid=ELLIPSOID)
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
    N = a_axis / np.sqrt(1.0 - ecc1**2.0 * np.cos(gridth) ** 2.0)

    # Coefficient for calculating Stokes coefficients from pressure field
    # SH Degree dependent factors with indirect loading components
    factors = gravtk.units(lmax=LMAX, a_axis=100.0 * a_axis, flat=flat)
    dfactor = factors.spatial(*LOVE).mmwe
    # use an integration factor for gridded data or
    # calculate from sin(theta)*dtheta*dphi
    int_fact = np.zeros((nlat))
    if WEIGHT is not None:
        # Weighting function for integrating gridded data
        int_fact[:] = np.broadcast_to(np.atleast_1d(WEIGHT), nlat)
    else:
        # Multiplying sin(th) with differentials of theta and phi
        # to calculate the integration factor at each latitude
        dphi = np.abs(phi[1] - phi[0])
        dth = np.abs(th[1] - th[0])
        int_fact[:] = np.sin(th) * dphi * dth

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    mm = np.arange(MMAX + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->mp...', mm, phi))

    # added option to precompute plms to improve computational speed
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # Fully-normalized Legendre Polynomials
    # Multiplying by integration factors [sin(theta)*dtheta*dphi]
    plm = np.einsum(
        'lmh...,h...->lmh...', PLM[: LMAX + 1, : MMAX + 1, :], int_fact
    )

    # gravitational acceleration at the Earth's mean spherical surface
    g0 = 9.80665
    # gravitational acceleration at the equator and at mean sea level
    ge = 9.780356
    # gravitational acceleration at the mean sea level over gridth
    gs = ge * (
        1.0
        + 5.2885e-3 * np.cos(gridth) ** 2
        - 5.9e-6 * np.cos(2.0 * gridth) ** 2
    )
    # calculate radii and gravity for each pressure level
    R = np.zeros((nlevels, nlat, nlon))
    gamma_h = np.zeros((nlevels, nlat, nlon))
    for p in range(nlevels):
        # orthometric height from List (1958)
        # as described in Boy and Chao (2005)
        orthometric = (1.0 - 0.002644 * np.cos(2.0 * gridth)) * GPH[p, :, :] + (
            1.0 - 0.0089 * np.cos(2.0 * gridth)
        ) * (GPH[p, :, :] ** 2) / 6.245e6
        # calculate X, Y and Z from geodetic latitude and longitude
        X = (N + GEOID + orthometric) * np.sin(gridth) * np.cos(gridphi)
        Y = (N + GEOID + orthometric) * np.sin(gridth) * np.sin(gridphi)
        Z = (N * (1.0 - ecc1**2.0) + GEOID + orthometric) * np.cos(gridth)
        # calculate radius of level
        R[p, :, :] = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
        # calculate normal gravity at each height above mean sea level
        gamma_h[p, :, :] = gs * (
            1.0
            - 2.0
            * (1.006803 - 0.06706 * np.cos(gridth) ** 2)
            * (orthometric / rad_e)
            + 3.0 * (orthometric / rad_e) ** 2
        )

    # total pressure factor
    pfactor = np.empty((nlat, nlon))
    # Initializing output spherical harmonic matrices
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX + 1, MMAX + 1))
    Ylms.slm = np.zeros((LMAX + 1, MMAX + 1))
    for l in range(0, LMAX + 1):  # equivalent to 0:LMAX
        mm = np.min([MMAX, l])  # truncate to MMAX if specified (if l > MMAX)
        m = slice(0, mm + 1)  # mm+1 elements between 0 and mm
        # zero out pressure factor for degree l
        pfactor[:, :] = 0.0
        # iterate over pressure levels
        for p in range(nlevels):
            # if using Swenson and Wahr (2002) or Boy and Chao (2005)
            if METHOD == 'SW02':
                # calculate pressure change/gravity ratio
                PG = pressure[p, :, :] / g0
                # add to pressure factor (pfactor) to integrate over levels
                pfactor += PG * np.power(
                    rad_e / (rad_e - GPH[p, :, :]) + (GEOID / rad_e), (l + 4)
                )
            elif METHOD == 'BC05':
                # calculate pressure change/gravity ratio
                PG = pressure[p, :, :] / gamma_h[p, :, :]
                # add to pressure factor (pfactor) to integrate over levels
                pfactor += PG * np.power(R[p, :, :] / rad_e, (l + 2))
        # Multiplying gridded data with sin/cos of m#phis
        # This will sum through all phis in the dot product
        # need to reform pfactor to lonXlat as is originally latXlon
        # output [m,theta]
        d = np.einsum('mp...,hp...->mh...', m_phi, -pfactor)
        # Summing product of plms and data over all latitudes
        ylm = np.einsum('mh...,mh...->m...', plm[l, m, :], d[m, :])
        # Multiplying by coefficients to normalize
        Ylms.clm[l, m] = dfactor[l] * ylm.real
        Ylms.slm[l, m] = dfactor[l] * ylm.imag

    # return the harmonics object
    return Ylms
