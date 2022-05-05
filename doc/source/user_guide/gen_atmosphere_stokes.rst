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

.. autofunction:: model_harmonics.gen_atmosphere_stokes
