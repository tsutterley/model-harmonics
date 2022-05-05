=====================
gldas_read_monthly.py
=====================

- Reads GLDAS monthly datafiles to calculate anomalies in total water storage from soil moisture, snow water equivalent and total canopy storage [Rodell2004]_ [Syed2008]_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_read_monthly.py

Calling Sequence
################

.. argparse::
    :filename: ../../GLDAS/gldas_read_monthly.py
    :func: arguments
    :prog: gldas_read_monthly.py
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

References
##########

.. [Rodell2004] M. Rodell et al., "The Global Land Data Assimilation System", *Bulletin of the American Meteorological Society*, 85(3), 381--394, (2004). `doi: 10.1175/BAMS-85-3-381 <https://doi.org/10.1175/BAMS-85-3-381>`_

.. [Syed2008] T. H. Syed et al., "Analysis of terrestrial water storage changes from GRACE and GLDAS", *Water Resources Research*, 44(2), (2008). `doi: 10.1029/2006WR005779 <https://doi.org/10.1029/2006WR005779>`_
