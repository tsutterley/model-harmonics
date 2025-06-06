#!/usr/bin/env python
u"""
gen_pressure_stokes.py
Written by Tyler Sutterley (03/2023)
Calculates spherical harmonic fields from spatial pressure fields

CALLING SEQUENCE:
    Ylms = gen_pressure_stokes(P, G, R, lon, lat, LMAX=60,
        PLM=PLM, LOVE=(hl,kl,ll))

INPUTS:
    P: Pressure [Pa]
    G: Gravitational acceleration [m/s^2]
    R: Radius at point [m]
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
    PLM: input Legendre polynomials
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    associated_legendre.py: Computes fully normalized associated
        Legendre polynomials
    units.py: class for converting spherical harmonic data to specific units
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
    Updated 03/2023: improve typing for variables in docstrings
    Updated 01/2023: refactored associated legendre polynomials
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 05/2021: define int/float precision to prevent deprecation warning
    Updated 02/2021: separate pressure and gravitational acceleration inputs
    Updated 01/2021: use harmonics class for spherical harmonic operations
    Updated 07/2020: added function docstrings
    Updated 04/2020: made Legendre polynomials and Love numbers options
        using the units class for converting to normalized spherical harmonics
    Updated 10/2018: separated into a single function for use with the
        ocean bottom pressure/atmospheric reanalysis/geocenter programs
    Updated 03/2018: simplified love number extrapolation if LMAX > 696
    Written 03/2018
"""

import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.associated_legendre import plm_holmes

# PURPOSE: calculates spherical harmonic fields from pressure fields
def gen_pressure_stokes(P, G, R, lon, lat, LMAX=60, MMAX=None,
    PLM=None, LOVE=None):
    r"""
    Converts pressure fields from the spatial domain to spherical
    harmonic coefficients :cite:p:`Boy:2005el,Swenson:2002kf`

    Parameters
    ----------
    P: np.ndarray
        Pressure (Pa)
    G: np.ndarray
        Gravitational acceleration (m/s\ :sup:`2`)
    R: np.ndarray
        Radius at point (m)
    lon: np.ndarray
        longitude array
    lat: np.ndarray
        latitude array
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    PLM: np.ndarray or NoneType, default None
        Legendre polynomials
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

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

    # grid dimensions
    nlat = np.int64(len(lat))
    # grid step
    dlon = np.abs(lon[1]-lon[0])
    dlat = np.abs(lat[1]-lat[0])
    # longitude degree spacing in radians
    dphi = dlon*np.pi/180.0
    # colatitude degree spacing in radians
    dth = dlat*np.pi/180.0

    # reformatting longitudes to range 0:360 (if previously -180:180)
    lon = np.squeeze(lon.copy())
    if np.any(lon < 0):
        lon_ind, = np.nonzero(lon < 0)
        lon[lon_ind] += 360.0
    # Longitude in radians
    phi = lon[np.newaxis,:]*np.pi/180.0
    # Colatitude in radians
    th = (90.0 - np.squeeze(lat.copy()))*np.pi/180.0

    # For gridded data: dmat = original data matrix
    sz = np.shape(P)
    # reforming data to lonXlat if input latXlon
    if (sz[0] == nlat):
        P = np.transpose(P)
        G = np.transpose(G)
        R = np.transpose(R)

    # Coefficient for calculating Stokes coefficients from pressure field
    # extract arrays of kl, hl, and ll Love Numbers
    factors = gravity_toolkit.units(lmax=LMAX).spatial(*LOVE)
    # Earth Parameters
    # Average Radius of the Earth [m]
    rad_e = factors.rad_e/100.0
    # SH Degree dependent factors with indirect loading components
    dfactor = factors.mmwe

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    m = np.arange(MMAX+1)
    ccos = np.cos(np.dot(m[:,np.newaxis],phi))
    ssin = np.sin(np.dot(m[:,np.newaxis],phi))

    # Calculates fully-normalized Legendre Polynomials with plm_holmes.py
    # Output is plm[l,m,th]
    plm = np.zeros((LMAX+1,MMAX+1,nlat))
    # added option to precompute plms to improve computational speed
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = plm_holmes(LMAX, np.cos(th))

    # Multiplying by integration factors [sin(theta)*dtheta*dphi]
    # truncate legendre polynomials to spherical harmonic order MMAX
    m = np.arange(MMAX+1)
    for j in range(0,nlat):
        plm[:,m,j] = PLM[:,m,j]*np.sin(th[j])*dphi*dth

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
        # Multiplying gridded data with sin/cos of m#phis
        # This will sum through all phis in the dot product
        # output [m,theta]
        pfactor = (P/G)*(R/rad_e)**(l+2)
        dcos = np.dot(ccos,pfactor)
        dsin = np.dot(ssin,pfactor)
        # Summing product of plms and data over all latitudes
        # axis=1 signifies the direction of the summation (colatitude (th))
        # ycos and ysin are the SH coefficients before normalizing
        yclm[l,m] = np.sum(plm[l,m,:]*dcos[m,:], axis=1)
        yslm[l,m] = np.sum(plm[l,m,:]*dsin[m,:], axis=1)
        # Multiplying by factors to normalize
        Ylms.clm[l,m] = dfactor[l]*yclm[l,m]
        Ylms.slm[l,m] = dfactor[l]*yslm[l,m]

    # return the harmonics object
    return Ylms
