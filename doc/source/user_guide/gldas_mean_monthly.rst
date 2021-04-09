gldas_mean_monthly.py
=====================

- Reads GLDAS monthly datafiles to calculate the multi-annual mean total water storage from soil moisture, snow water equivalent and total canopy storage [Rodell2004]_ [Syed2008]_

Calling Sequence
################

.. code-block:: bash

    python gldas_mean_monthly.py --directory <path_to_directory> --version 2.1 CLSM NOAH VIC

`Source code`__

.. __ : https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_mean_monthly.py

Inputs
######

- ``'CLM'``: GLDAS Common Land Model
- ``'CLSM'``: GLDAS Catchment Land Surface Model
- ``'MOS'``: GLDAS Mosaic model
- ``'NOAH'``: GLDAS Noah model
- ``'VIC'``: GLDAS Variable Infiltration Capacity model

Command Line Options
####################

- ``-D X``, ``--directory X``: Working data directory
- ``-M X``, ``--mean X``: Year range for mean
- ``-S X``, ``--spacing X``: Spatial resolution of models to run

    * ``'10'``: 1.0 degrees latitude/longitude
    * ``'025'``: 0.25 degrees latitude/longitude
- ``-v X``, ``--version X``: GLDAS model version to run
- ``-F X``, ``--format X``: input and output data format

    * ``'ascii'``
    * ``'netcdf'``
    * ``'HDF5'``
- ``-M X``, ``--mode X``: Permission mode of directories and files
- ``-V``, ``--verbose``: Output information for each output file

References
##########

.. [Rodell2004] M. Rodell et al., "The Global Land Data Assimilation System", *Bulletin of the American Meteorological Society*, 85(3), 381--394, (2004). `doi: 10.1175/BAMS-85-3-381 <https://doi.org/10.1175/BAMS-85-3-381>`_

.. [Syed2008] T. H. Syed et al., "Analysis of terrestrial water storage changes from GRACE and GLDAS", *Water Resources Research*, 44(2), (2008). `doi: 10.1029/2006WR005779 <https://doi.org/10.1029/2006WR005779>`_
