======================
gen_pressure_stokes.py
======================

- Converts pressure fields from the spatial domain to spherical harmonic coefficients [Boy2005]_ [Swenson2002]_

Calling Sequence
################

.. code-block:: python

    from model_harmonics.gen_pressure_stokes import gen_pressure_stokes
    from gravity_toolkit.plm_holmes import plm_holmes
    PLM, dPLM = plm_holmes(LMAX, np.cos(th))
    Ylms = gen_pressure_stokes(P, G, R, lon, lat, LMAX=LMAX, PLM=PLM, LOVE=(hl,kl,ll))

`Source code`__

.. __: https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/gen_pressure_stokes.py

.. autofunction:: model_harmonics.gen_pressure_stokes
