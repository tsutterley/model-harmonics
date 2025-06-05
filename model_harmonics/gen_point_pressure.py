#!/usr/bin/env python
u"""
gen_point_pressure.py
Written by Tyler Sutterley (03/2023)
Calculates gravitational spherical harmonic coefficients for pressure
    values at individual points assuming a disc geometry

CALLING SEQUENCE:
    Ylms = gen_point_pressure(P, G, R, lon, lat, LMAX=LMAX)

INPUTS:
    P: Pressure [Pa]
    G: Gravitational acceleration [m/s^2]
    R: Radius at point [m]
    lon: longitude of points
    lat: latitude of points

OUTPUTS:
    clm: cosine spherical harmonic coefficients (geodesy normalization)
    slm: sine spherical harmonic coefficients (geodesy normalization)
    l: spherical harmonic degree to LMAX
    m: spherical harmonic order to MMAX

OPTIONS:
    AREA: Area of each pressure cell [m^2]
    LMAX: Upper bound of Spherical Harmonic Degrees
    MMAX: Upper bound of Spherical Harmonic Orders
    LOVE: input load Love numbers up to degree LMAX (hl,kl,ll)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    scipy: Scientific Tools for Python (https://docs.scipy.org/doc/)

PROGRAM DEPENDENCIES:
    legendre.py: Computes associated Legendre polynomials for degree l
    units.py: class for converting spherical harmonic data to specific units
    harmonics.py: spherical harmonic data class for processing GRACE/GRACE-FO
    destripe_harmonics.py: calculates the decorrelation (destriping) filter
        and filters the GRACE/GRACE-FO coefficients for striping errors

REFERENCES:
    I. M. Longman, Journal of Geophysical Research, 67(2), 1962
        https://doi.org/10.1029/JZ067i002p00845
    W. E. Farrell, Reviews of Geophysics and Space Physics, 10(3), 1972
        https://doi.org/10.1029/RG010i003p00761
    H. N. Pollack, Journal of Geophysical Research, 78(11), 1973
        https://doi.org/10.1029/JB078i011p01760
    T. Jacob et al., Journal of Geodesy, 86, 337-358, 2012
        https://doi.org/10.1007/s00190-011-0522-7

UPDATE HISTORY:
    Updated 03/2023: simplified recursion and unit degree factors
        improve typing for variables in docstrings
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 02/2021
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.legendre import legendre

def gen_point_pressure(P, G, R, lon, lat, AREA=None, LMAX=60, MMAX=None,
    LOVE=None):
    r"""
    Calculates gravitational spherical harmonic coefficients for pressure
    values at individual points assuming a disc geometry
    :cite:p:`Boy:2005el,Longman:1962ev,Farrell:1972cm,Pollack:1973gi,Swenson:2002kf`


    Parameters
    ----------
    P: np.ndarray
        Pressure (Pa)
    G: np.ndarray
        Gravitational acceleration (m/s\ :sup:`2`)
    R: np.ndarray
        Radius at point (m)
    lon: np.ndarray
        longitude of points
    lat: np.ndarray
        latitude of points
    AREA: np.ndarray or NoneType, default None
        Area of each pressure cell (m\ :sup:`2`)
    LMAX: int, default 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
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

    # upper bound of spherical harmonic orders (default == LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # convert output longitude and latitude into radians
    npts = len(lon.flatten())
    phi = np.pi*lon.flatten()/180.0
    theta = np.pi*(90.0 - lat.flatten())/180.0

    # SH Degree dependent factors to convert into fully normalized SH's
    factors = gravity_toolkit.units(lmax=LMAX).spatial(*LOVE)
    # Earth Parameters
    # Average Radius of the Earth [m]
    rad_e = factors.rad_e/100.0
    # Coefficient for calculating Stokes coefficients for a disc load
    # From Jacob et al (2012), Farrell (1972) and Longman (1962)
    dfactor = 4.0*np.pi*factors.mmwe/(1.0 + 2.0*factors.l)

    # Calculating legendre polynomials of the disc
    # alpha will be 1 - the ratio of the input area with the half sphere
    alpha = (1.0 - AREA.flatten()/(2.0*np.pi*rad_e**2))
    # seeds for Legendre Polynomial recursion (degrees l-1, l)
    Pm1 = np.ones((npts))
    Pl = np.ones((npts))

    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    # for each degree l
    for l in range(LMAX+1):
        m1 = np.min([l,MMAX]) + 1
        # Calculating legendre polynomials for degree l+1
        Pp1 = ((2.0*l + 1.0)/(l+1.0))*alpha*Pl - (l/(l + 1.0))*Pm1
        # legendre polynomials of the disc (unnormalized)
        # from Longman (1962) and Jacob et al (2012)
        Pdisc = (Pm1 - Pp1)/2.0
        # calculate pressure/gravity ratio for all points
        # convolve with legendre polynomials of the disc
        # and the radius ratios
        PGR = Pdisc*(P.flatten()/G.flatten())*(R.flatten()/rad_e)**(l+2)
        SPH = spherical_harmonic_matrix(l,PGR,phi,theta,dfactor[l])
        # truncate to spherical harmonic order and save to output
        Ylms.clm[l,:m1] = SPH.real[:m1]
        Ylms.slm[l,:m1] = SPH.imag[:m1]
        # update unnormalized Legendre polynomials for recursion
        Pm1[:] = np.copy(Pl)
        Pl[:] = np.copy(Pp1)
    # return the output spherical harmonics object
    return Ylms

# calculate spherical harmonics of degree l evaluated at (theta,phi)
def spherical_harmonic_matrix(l, PGR, phi, theta, coeff):
    """
    Calculates spherical harmonics of degree l evaluated at coordinates

    Parameters
    ----------
    l: int
        spherical harmonic degree
    PGR: np.ndarray
        pressure/gravity ratio
    phi: np.ndarray
        longitude of points in radians
    theta: np.ndarray
        colatitude of points in radians
    coeff: np.ndarray
        degree-dependent factor for converting units

    Returns
    -------
    Ylms: np.ndarray
        spherical harmonic coefficients in Eulerian form
    """
    # calculate normalized legendre polynomials (points, order)
    Pl = legendre(l, np.cos(theta), NORMALIZE=True).T
    # spherical harmonic orders up to degree l
    m = np.arange(0,l+1)
    # calculate Euler's of spherical harmonic order multiplied by azimuth phi
    mphi = np.exp(1j*np.dot(np.squeeze(phi)[:,np.newaxis],m[np.newaxis,:]))
    # reshape pressure/gravity ratio to order
    D = np.kron(np.ones((1,l+1)), PGR[:,np.newaxis])
    # calculate spherical harmonics and multiply by coefficients and data
    Ylms = coeff*D*Pl*mphi
    # calculate the sum over all points and return harmonics for degree l
    return np.sum(Ylms,axis=0)
