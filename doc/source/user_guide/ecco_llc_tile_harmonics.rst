==========================
ecco_llc_tile_harmonics.py
==========================

- Reads monthly ECCO ocean bottom pressure anomalies from LLC tiles and converts to spherical harmonic coefficients [Boy2005]_

Calling Sequence
################

.. code-block:: bash

    python ecco_llc_tile_harmonics.py --directory <path_to_directory> V5alpha

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_llc_tile_harmonics.py

Inputs
######

- ECCO Version 4 or 5 models

    * ``'V4r4'``: Version 4, Revision 4
    * ``'V5alpha'``: ECCO Version 5, Alpha release

Command Line Options
####################

- ``-D X``, ``--directory X``: working data directory
- ``-Y X``, ``--year X``: Years to run
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
- ``-F X``, ``--format X``: output data format

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
- ``-V``, ``--verbose``: verbose output of processing run
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Boy2005] J.-P. Boy and B. F. Chao, "Precise evaluation of atmospheric loading effects on Earth's time‐variable gravity field", *Journal of Geophysical Research: Solid Earth*, 110(B08412), (2005). `doi: 10.1029/2002JB002333 <https://doi.org/10.1029/2002JB002333>`_

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
