===========================
model_level_coefficients.py
===========================

- Creates a netCDF4 file of reanalysis A and B coefficients for model levels
- Model level coefficients are obtained using equation 3.17 of [Simmons1981]_ and the methodology of [Trenberth1993]_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/model_level_coefficients.py

Calling Sequence
################

.. argparse::
    :filename: ../../reanalysis/model_level_coefficients.py
    :func: arguments
    :prog: model_level_coefficients.py
    :nodescription:
    :nodefault:

    model : @after
        * `ERA-Interim <http://apps.ecmwf.int/datasets/data/interim-full-moda>`_
        * `ERA5 <http://apps.ecmwf.int/data-catalogues/era5/?class=ea>`_
        * `MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_

References
##########

.. [Simmons1981] A. J. Simmons and D. M. Burridge, "An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates" *Monthly Weather Review*, 109(4), 758--766, (1981). `doi: 10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2`__

.. __: https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2

.. [Trenberth1993] K. E. Trenberth, J. C. Berry, and L. E. Buja, "Vertical Interpolation and Truncation of Model-coordinate Data", Technical Note, No. NCAR/TN-396+STR, (1993). `doi: 10.5065/D6HX19NH <https://doi.org/10.5065/D6HX19NH>`_
