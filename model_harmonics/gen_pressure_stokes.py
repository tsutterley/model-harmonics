#!/usr/bin/env python
"""
gen_pressure_stokes.py
Written by Tyler Sutterley (07/2026)
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
    Updated 07/2026: use np.einsum for spherical harmonic summations
        use np.radians to convert from degrees to radians
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
import gravity_toolkit as gravtk


# PURPOSE: calculates spherical harmonic fields from pressure fields
def gen_pressure_stokes(
    P, G, R, lon, lat, LMAX=60, MMAX=None, PLM=None, LOVE=None
):
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

    # calculate longitude and colatitude arrays in radians
    phi = np.radians(np.squeeze(lon))
    th = np.radians(90.0 - np.squeeze(lat))
    # reformatting longitudes to range 0:360 (if previously -180:180)
    phi = np.where(phi < 0, phi + 2.0 * np.pi, phi)
    # grid step in radians
    dphi = np.abs(phi[1] - phi[0])
    dth = np.abs(th[1] - th[0])

    # For gridded data: dmat = original data matrix
    sz = np.shape(P)
    # reforming data to lonXlat if input latXlon
    if sz[0] == len(lat):
        P = np.transpose(P)
        G = np.transpose(G)
        R = np.transpose(R)

    # Coefficient for calculating Stokes coefficients from pressure field
    # extract arrays of kl, hl, and ll Love Numbers
    factors = gravtk.units(lmax=LMAX).spatial(*LOVE)
    # Earth Parameters
    # Average Radius of the Earth [m]
    rad_e = factors.rad_e / 100.0
    # SH Degree dependent factors with indirect loading components
    dfactor = factors.mmwe
    # Multiplying sin(th) with differentials of theta and phi
    # to calculate the integration factor at each latitude
    int_fact = np.sin(th) * dphi * dth

    # Calculating cos/sin of phi arrays
    # output [m,phi]
    mm = np.arange(MMAX + 1)
    m_phi = np.exp(1j * np.einsum('m...,p...->mp...', mm, phi))

    # Calculate polynomials using Holmes and Featherstone (2002) relation
    if PLM is None:
        # if plms are not pre-computed: calculate Legendre polynomials
        PLM, dPLM = gravtk.plm_holmes(LMAX, np.cos(th))

    # Fully-normalized Legendre Polynomials
    # Multiplying by integration factors [sin(theta)*dtheta*dphi]
    plm = np.einsum(
        'lmh...,h...->lmh...', PLM[: LMAX + 1, : MMAX + 1, :], int_fact
    )

    # Initializing output spherical harmonic matrices
    Ylms = gravtk.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX + 1, MMAX + 1))
    Ylms.slm = np.zeros((LMAX + 1, MMAX + 1))
    for l in range(0, LMAX + 1):  # equivalent to 0:LMAX
        mm = np.min([MMAX, l])  # truncate to MMAX (if l > MMAX)
        m = slice(0, mm + 1)  # mm+1 elements between 0 and mm
        # Multiplying gridded data with sin/cos of m#phis
        # This will sum through all phis in the dot product
        # output [m,theta]
        pfactor = (P / G) * np.power(R / rad_e, (l + 2))
        d = np.einsum('mp...,ph...->mh...', m_phi, pfactor)
        # Summing product of plms and data over all latitudes
        ylm = np.einsum('mh...,mh...->m...', plm[l, m, :], d[m, :])
        # Multiplying by factors to normalize
        Ylms.clm[l, m] = dfactor[l] * ylm.real
        Ylms.slm[l, m] = dfactor[l] * ylm.imag

    # return the harmonics object
    return Ylms
