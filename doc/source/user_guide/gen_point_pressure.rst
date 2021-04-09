=====================
gen_point_pressure.py
=====================

- Calculates gravitational spherical harmonic coefficients for pressure values at individual points assuming a disc geometry [Boy2005]_ [Swenson2002]_

Calling Sequence
################

.. code-block:: python

    from model_harmonics.gen_point_pressure import gen_point_pressure
    Ylms = gen_point_pressure(P, G, R, lon, lat, AREA=AREA, LMAX=LMAX, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/model_harmonics/gen_point_pressure.py

Arguments
#########

1. ``P``: Pressure (Pa)
2. ``G``: Gravitational acceleration (m/s\ :sup:`2`)
3. ``R``: Earth's radius at each data point (m)
4. ``lon``: longitude of points
5. ``lat``: latitude of points

Keyword arguments
#################

- ``AREA``: Area of each pressure cell (m\ :sup:`2`)
- ``LMAX``:  maximum spherical harmonic degree of the output harmonics
- ``MMAX``: maximum spherical harmonic order of the output harmonics
- ``LOVE``: input load Love numbers up to degree of truncation

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
