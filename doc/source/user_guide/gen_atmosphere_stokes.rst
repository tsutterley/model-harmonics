========================
gen_atmosphere_stokes.py
========================

 - Converts 3D atmospheric geopotential height and pressure difference fields from the spatial domain to spherical harmonic coefficients [Boy2005]_ [Swenson2002]_

 Calling Sequence
 ################

 .. code-block:: python

    from model_harmonics.gen_atmosphere_stokes import gen_atmosphere_stokes
    from gravity_toolkit.plm_holmes import plm_holmes
    PLM,dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_atmosphere_stokes(GPH, pressure, lon, lat, LMAX=LMAX,
        ELLIPSOID='WGS84', GEOID=geoid, PLM=PLM, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/model_harmonics/gen_atmosphere_stokes.py

Arguments
#########

1. ``GPH``: geopotential heights at model levels
2. ``pressure``: pressure differences between model levels
3. ``lon``: longitude array
4. ``lat``: latitude array

Keyword arguments
#################

- ``LMAX``:  maximum spherical harmonic degree of the output harmonics
- ``MMAX``: maximum spherical harmonic order of the output harmonics
- ``ELLIPSOID``: reference ellipsoid name
- ``GEOID``: geoid height
- ``PLM``: input Legendre polynomials
- ``LOVE``: input load Love numbers up to degree of truncation
- ``METHOD``: method of integrating over pressure levels

    * ``'SW02'``: Swenson and Wahr (2002) [Swenson2002]_
    * ``'BC05'``: Boy and Chao (2005) [Boy2005]_

Returns
#######

- ``clm``: Cosine spherical harmonic coefficients (geodesy normalization)
- ``slm``: Sine spherical harmonic coefficients (geodesy normalization)
- ``l``: spherical harmonic degree to ``LMAX``
- ``m``: spherical harmonic order to ``MMAX``

References
##########

.. [Boy2005] J.-P. Boy and B. F. Chao, "Precise evaluation of atmospheric loading effects on Earth's time‐variable gravity field", *Journal of Geophysical Research: Solid Earth*, 110(B08412), (2005). `doi: 10.1029/2002JB002333 <https://doi.org/10.1029/2002JB002333>`_

.. [Swenson2002] S. Swenson and J. Wahr, "Estimated effects of the vertical structure of atmospheric mass on the time‐variable geoid", *Journal of Geophysical Research*, 107(B9), 2194, (2002) `doi: 10.1029/2000JB000024 <https://doi.org/10.1029/2000JB000024>`_
