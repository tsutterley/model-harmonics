#!/usr/bin/env python
u"""
gen_point_pressure.py
Written by Tyler Sutterley (04/2022)
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
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 02/2021
"""
import numpy as np
import gravity_toolkit.units
import gravity_toolkit.harmonics
from gravity_toolkit.legendre import legendre

def gen_point_pressure(P, G, R, lon, lat, AREA=None, LMAX=60, MMAX=None,
    LOVE=None):
    """
    Calculates gravitational spherical harmonic coefficients for pressure
        values at individual points assuming a disc geometry

    Parameters
    ----------
    P: float
        Pressure (Pa)
    G: float
        Gravitational acceleration (m/s\ :sup:`2`)
    R: float
        Radius at point (m)
    lon: float
        longitude of points
    lat: float
        latitude of points
    AREA: float or NoneType, default None
        Area of each pressure cell (m\ :sup:`2`)
    LMAX: int, defualt 60
        Upper bound of Spherical Harmonic Degrees
    MMAX: int or NoneType, default None
        Upper bound of Spherical Harmonic Orders
    LOVE: tuple or NoneType, default None
        Load Love numbers up to degree LMAX (``hl``, ``kl``, ``ll``)

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
    .. [Farrell1972] W. E. Farrell, "Deformation of the Earth by surface loads",
        *Reviews of Geophysics and Space Physics*, 10(3), (1972).
        `doi: 10.1029/RG010i003p00761 <https://doi.org/10.1029/RG010i003p00761>`_
    .. [Jacob2012] T. Jacob et al., "Estimating geoid height change in North America:
        past, present and future", *Journal of Geodesy*, 86, 337-358, (2012).
        `doi: 10.1007/s00190-011-0522-7 <https://doi.org/10.1007/s00190-011-0522-7>`_
    .. [Longman1962] I. M. Longman, "A Green's function for determining
        the deformation of the Earth under surface mass loads: 1. Theory",
        *Journal of Geophysical Research*, 67(2), (1962).
        `doi: 10.1029/JZ067i002p00845 <https://doi.org/10.1029/JZ067i002p00845>`_
    .. [Pollack1973] H. N. Pollack, "Spherical harmonic representation of the
        gravitational potential of a point mass, a spherical cap, and a
        spherical rectangle", *Journal of Geophysical Research*, 78(11), (1973).
        `doi: 10.1029/JB078i011p01760 <https://doi.org/10.1029/JB078i011p01760>`_
    .. [Swenson2002] S. Swenson and J. Wahr, "Estimated effects of the vertical
        structure of atmospheric mass on the time‐variable geoid",
        *Journal of Geophysical Research*, 107(B9), 2194, (2002).
        `doi: 10.1029/2000JB000024 <https://doi.org/10.1029/2000JB000024>`_
    """

    # upper bound of spherical harmonic orders (default == LMAX)
    if MMAX is None:
        MMAX = np.copy(LMAX)

    # convert output longitude and latitude into radians
    npts = len(lon.flatten())
    phi = np.pi*lon.flatten()/180.0
    theta = np.pi*(90.0 - lat.flatten())/180.0

    # extract arrays of kl, hl, and ll Love Numbers
    hl,kl,ll = LOVE
    # SH Degree dependent factors to convert into fully normalized SH's
    factors = gravity_toolkit.units(lmax=LMAX).spatial(hl,kl,ll)
    # Earth Parameters
    # Average Density of the Earth [kg/m^3]
    rho_e = 1000.0*factors.rho_e
    # Average Radius of the Earth [m]
    rad_e = factors.rad_e/100.0
    # Coefficient for calculating Stokes coefficients for a disc load
    # From Jacob et al (2012), Farrell (1972) and Longman (1962)
    dfactor = 3.0*(1.0 + kl[factors.l])/(rad_e*rho_e*(1.0 + 2.0*factors.l)**2)

    # Calculating legendre polynomials of the disc
    # alpha will be 1 - the ratio of the input area with the half sphere
    alpha = (1.0 - AREA.flatten()/(2.0*np.pi*rad_e**2))
    # seed for Legendre Polynomial recursion
    Pm2 = np.copy(alpha)
    Pm1 = np.ones((npts))

    # Initializing output spherical harmonic matrices
    Ylms = gravity_toolkit.harmonics(lmax=LMAX, mmax=MMAX)
    Ylms.clm = np.zeros((LMAX+1,MMAX+1))
    Ylms.slm = np.zeros((LMAX+1,MMAX+1))
    # for each degree l
    for l in range(LMAX+1):
        m1 = np.min([l,MMAX]) + 1
        # unnormalized Legendre polynomials for degree l
        if (l == 0):
            # l=0 is a special case
            Pl = np.ones((npts))
        else:
            # Calculating legendre polynomials for degree l
            Pl = ((2.0*l-1.0)/l)*alpha*Pm1 - ((l-1.0)/l)*Pm2
        # Calculating legendre polynomials for degree l+1
        Pp1 = ((2.0*l+1.0)/(l+1.0))*alpha*Pl - (l/(l+1.0))*Pm1
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
        Pm2[:] = np.copy(Pm1)
        Pm1[:] = np.copy(Pl)
    # return the output spherical harmonics object
    return Ylms

# calculate spherical harmonics of degree l evaluated at (theta,phi)
def spherical_harmonic_matrix(l,PGR,phi,theta,coeff):
    """
    Calculates spherical harmonics of degree l evaluated at coordinates

    Parameters
    ----------
    l: int
        spherical harmonic degree
    PGR: float
        pressure/gravity ratio
    phi: float
        longitude of points in radians
    theta: float
        colatitude of points in radians
    coeff: float
        degree-dependent factor for converting units

    Returns
    -------
    Ylms: float
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
