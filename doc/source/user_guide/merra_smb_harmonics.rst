======================
merra_smb_harmonics.py
======================

- Reads monthly MERRA-2 surface mass balance anomalies and converts to spherical harmonic coefficients [Wahr1998]_

Calling Sequence
################

.. code-block:: bash

    python merra_smb_harmonics.py --directory <path_to_directory> SMB

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_harmonics.py

Inputs
######

- ``'SMB'``: Surface Mass Balance
- ``'ACCUM'``: Snowfall accumulation
- ``'PRECIP'``: Total Precipitation
- ``'RAINFALL'``: Total Rainfall
- ``'SUBLIM'``: Evaporation and Sublimation
- ``'RUNOFF'``: Meltwater Runoff

Command Line Options
####################

- ``-D X``, ``--directory X``: working data directory
- ``-m X``, ``--mean X``: Year range for mean
- ``-Y X``, ``--year X``: Years to run
- ``-R X``, ``--region X``: region name for subdirectory
- ``--mask X``: netCDF4 mask file for reducing to regions
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``-n X``, ``--love X``: Load Love numbers dataset

    * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
    * ``1``: Gegout (2005) values from PREM [Gegout2010]_
    * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``--reference X``: Reference frame for load love numbers

    * ``'CF'``: Center of Surface Figure (default)
    * ``'CM'``: Center of Mass of Earth System
    * ``'CE'``: Center of Mass of Solid Earth
- ``-F X``, ``--format X``: input and output data format

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
