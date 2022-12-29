#!/usr/bin/env python
u"""
constants.py
Written by Tyler Sutterley (12/2022)

Gravitational and ellipsoidal parameters

CALLING SEQUENCE
    wgs84 = constants('WGS84', units='MKS')

INPUT:
    ellipsoid - reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967 (IAU ellipsoid)
        GRS80 = Geodetic Reference System 1980
        WGS72 = World Geodetic System 1972
        WGS84 = World Geodetic System 1984
        ATS77 = Quasi-earth centred ellipsoid for ATS77
        NAD27 = North American Datum 1927
        NAD83 = North American Datum 1983
        INTER = International
        KRASS = Krassovsky (USSR)
        MAIRY = Modified Airy (Ireland 1965/1975)
        TOPEX = TOPEX/POSEIDON ellipsoid
        EGM96 = EGM 1996 gravity model
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Hofmann-Wellenhof and Moritz (2006)

UPDATE HISTORY:
    Forked 12/2022 from ref_ellipsoid.py
"""
import numpy as np

class constants(object):
    """
    Class for gravitational and ellipsoidal parameters

    Parameters
    ----------
    ellipsoid: str, default 'WGS84'
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model
    units: str, default `MKS`
        Output units

            - ``'MKS'``: meters, kilograms, seconds
            - ``'CGS'``: centimeters, grams, seconds

    References
    ----------
    .. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz,
        *Physical Geodesy*, 2nd Edition, 403 pp., (2006).
        `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
    """
    np.seterr(invalid='ignore')
    def __init__(self, ellipsoid='WGS84', units='MKS'):
        # validate units
        assert units in ('MKS', 'CGS')

        # set parameters for ellipsoid
        if ellipsoid.upper() in ('CLK66','NAD27'):
            # Clarke 1866
            self.a_axis = 6378206.4# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/294.9786982# flattening of the ellipsoid

        elif ellipsoid.upper() in ('GRS80','NAD83'):
            # Geodetic Reference System 1980
            # North American Datum 1983
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.26# flattening of the ellipsoid
            self.GM = 3.986005e14# [m^3/s^2] Geocentric Gravitational Constant

        elif (ellipsoid.upper() == 'GRS67'):
            # Geodetic Reference System 1967
            # International Astronomical Union (IAU ellipsoid)
            self.a_axis = 6378160.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.247167427# flattening of the ellipsoid
            self.GM = 3.98603e14# [m^3/s^2] Geocentric Gravitational Constant
            self.omega = 7292115.1467e-11# angular velocity of the Earth [rad/s]

        elif (ellipsoid.upper() == 'WGS72'):
            # World Geodetic System 1972
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.26# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'WGS84'):
            # World Geodetic System 1984
            self.a_axis = 6378137.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257223563# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'ATS77'):
            # Quasi-earth centred ellipsoid for ATS77
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'KRASS'):
            # Krassovsky (USSR)
            self.a_axis = 6378245.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.3# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'INTER'):
            # International
            self.a_axis = 6378388.0# [m] semimajor axis of the ellipsoid
            self.flat = 1/297.0# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'MAIRY'):
            # Modified Airy (Ireland 1965/1975)
            self.a_axis = 6377340.189# [m] semimajor axis of the ellipsoid
            self.flat = 1/299.3249646# flattening of the ellipsoid

        elif (ellipsoid.upper() == 'TOPEX'):
            # TOPEX/POSEIDON ellipsoid
            self.a_axis = 6378136.3# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257# flattening of the ellipsoid
            self.GM = 3.986004415e14# [m^3/s^2]

        elif (ellipsoid.upper() == 'EGM96'):
            # EGM 1996 gravity model
            self.a_axis = 6378136.3# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.256415099# flattening of the ellipsoid
            self.GM = 3.986004415e14# [m^3/s^2]

        elif (ellipsoid.upper() == 'HGH80'):
            # Hughes 1980 Ellipsoid used in some NSIDC data
            self.a_axis = 6378273.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.279411123064# flattening of the ellipsoid

        else:
            raise ValueError('Incorrect reference ellipsoid Name')

        # set default parameters if not listed as part of ellipsoid
        if ellipsoid.upper() not in ('GRS80','GRS67','NAD83','TOPEX','EGM96'):
            # for ellipsoids not listing the Geocentric Gravitational Constant
            self.GM = 3.986004418e14# [m^3/s^2]

        if ellipsoid.upper() not in ('GRS67'):
            # for ellipsoids not listing the angular velocity of the Earth
            self.omega = 7292115e-11# [rad/s]

        # universal gravitational constant [N*m^2/kg^2]
        self.G = 6.67430e-11

        # convert units to CGS
        if (units == 'CGS'):
            self.a_axis *= 100.0
            self.GM *= 1e6
            self.G *= 1000.0 # [dyn*cm^2/g^2]

    # mean radius of the Earth having the same volume
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 - f)
    @property
    def rad_e(self):
        """Average radius of the Earth with same volume as ellipsoid
        """
        return self.a_axis*(1.0 - self.flat)**(1.0/3.0)

    # semiminor axis of the ellipsoid
    @property
    def b_axis(self):
        """Semi-minor axis of the ellipsoid
        """
        return (1.0 - self.flat)*self.a_axis

    # Ratio between ellipsoidal axes
    @property
    def ratio(self):
        """Ratio between ellipsoidal axes
        """
        return (1.0 - self.flat)

    # Polar radius of curvature
    @property
    def rad_p(self):
        return self.a_axis/(1.0 - self.flat)

    # Linear eccentricity
    @property
    def ecc(self):
        """Linear eccentricity
        """
        return np.sqrt((2.0*self.flat - self.flat**2)*self.a_axis**2)

    # first numerical eccentricity
    @property
    def ecc1(self):
        """First numerical eccentricity
        """
        return self.ecc/self.a_axis

    # second numerical eccentricity
    @property
    def ecc2(self):
        """Second numerical eccentricity
        """
        return self.ecc/self.b_axis

    # m parameter [omega^2*a^2*b/(GM)]
    # p. 70, Eqn.(2-137)
    @property
    def m(self):
        """m Parameter
        """
        return self.omega**2*((1 - self.flat)*self.a_axis**3)/self.GM

    # q
    # p. 67, Eqn.(2-113)
    @property
    def q(self):
        """q Parameter
        """
        return 0.5*((1.0 + 3.0/(self.ecc2**2))*np.arctan(self.ecc2)-3.0/self.ecc2)

    # q_0
    # p. 67, Eqn.(2-113)
    @property
    def q0(self):
        """q\ :sub:`0`\ Parameter
        """
        return 3*(1.0 + 1.0/(self.ecc2**2)) * \
            (1.0 -1.0/self.ecc2*np.arctan(self.ecc2)) - 1.0

    # J_2 p. 75 Eqn.(2-167), p. 76 Eqn.(2-172)
    @property
    def J2(self):
        """Oblateness coefficient
        """
        return (self.ecc1**2)*(1.0 - 2.0*self.m*self.ecc2/(15.0*self.q))/3.0

    # Normalized C20 harmonic
    # p. 60, Eqn.(2-80)
    @property
    def C20(self):
        """Normalized C\ :sub:`20`\ harmonic
        """
        return -self.J2/np.sqrt(5.0)

    # Normal gravity at the equator
    # p. 71, Eqn.(2-141)
    @property
    def ga(self):
        """Normal gravity at the equator
        """
        return self.GM/(self.a_axis*self.b_axis) * \
            (1.0 - self.m - self.m*self.ecc2*self.q0/(6.0*self.q))

    # Normal gravity at the pole
    # p. 71, Eqn.(2-142)
    @property
    def gb(self):
        """Normal gravity at the pole
        """
        return self.GM/(self.a_axis**2.0) * \
            (1.0 + self.m*self.ecc2*self.q0/(3.0*self.q))

    # ratio between gravity at pole versus gravity at equator
    @property
    def dk(self):
        """Ratio between gravity at pole versus gravity at equator
        """
        return self.b_axis*self.gb/(self.a_axis*self.ga) - 1.0

    # Normal potential at the ellipsoid
    # p. 68, Eqn.(2-123)
    @property
    def U0(self):
        """Normal potential at the ellipsoid
        """
        return self.GM/self.ecc*np.arctan(self.ecc2) + \
            (1.0/3.0)*self.omega**2*self.a_axis**2

    # Surface area of the reference ellipsoid
    @property
    def area(self):
        """Surface area of the ellipsoid
        """
        return np.pi*self.a_axis**2.0 * \
            (2.0 + ((1.0 - self.ecc1**2)/self.ecc1) *
                np.log((1.0 + self.ecc1)/(1.0 - self.ecc1)))

    # Volume of the reference ellipsoid
    @property
    def volume(self):
        """Volume of the ellipsoid
        """
        return (4.0*np.pi/3.0)*(self.a_axis**3.0)*(1.0 - self.ecc1**2.0)**0.5

    # Average density
    @property
    def rho_e(self):
        """Average density
        """
        return self.GM/(self.G*self.volume)
