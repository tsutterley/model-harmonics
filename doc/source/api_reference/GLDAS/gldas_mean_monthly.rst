=====================
gldas_mean_monthly.py
=====================

- Reads GLDAS monthly datafiles to calculate the multi-annual mean total water storage from soil moisture, snow water equivalent and total canopy storage :cite:p:`Rodell:2004ke` :cite:p:`Syed:2008ia`

`Source code`__

.. __ : https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mean_monthly.py

Calling Sequence
################

.. argparse::
    :filename: gldas_mean_monthly.py
    :func: arguments
    :prog: gldas_mean_monthly.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'CLM'``: GLDAS Common Land Model
        * ``'CLSM'``: GLDAS Catchment Land Surface Model
        * ``'MOS'``: GLDAS Mosaic model
        * ``'NOAH'``: GLDAS Noah model
        * ``'VIC'``: GLDAS Variable Infiltration Capacity model

    --spacing -S : @after
        * ``'10'``: 1.0 degrees latitude/longitude
        * ``'025'``: 0.25 degrees latitude/longitude
