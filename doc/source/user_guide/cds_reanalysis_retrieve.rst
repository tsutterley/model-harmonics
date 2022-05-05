==========================
cds_reanalysis_retrieve.py
==========================

- Retrieves ERA5 reanalysis netCDF4 datasets from the CDS Web API

    * 2-metre Temperature (t2m)
    * Surface Pressure (ps)
    * Mean Sea Level Pressure (msl)
    * Temperature (t) and Specific Humidity (q) on Model Levels
    * Invariant Parameters

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/cds_reanalysis_retrieve.py

.. argparse::
    :filename: ../../reanalysis/cds_reanalysis_retrieve.py
    :func: arguments
    :prog: cds_reanalysis_retrieve.py
    :nodescription:
    :nodefault:

    --surface -S : @after
        * ``'MSL'``: mean sea level pressure field
        * ``'SP'``: surface pressure field
        * ``'T2m'``: 2-metre temperature field
        * ``'P-E'``: precipitation and evaporation fields
