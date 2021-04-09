===================================
reanalysis_atmospheric_harmonics.py
===================================

- Reads atmospheric geopotential heights fields from reanalysis and calculates sets of spherical harmonics using a 3D geometry [Boy2005]_ [Swenson2002]_

Calling Sequence
################

.. code-block:: bash

    python reanalysis_atmospheric_harmonics.py --directory <path_to_directory> ERA5 MERRA-2

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_atmospheric_harmonics.py

Inputs
######

- `ERA-Interim <http://apps.ecmwf.int/datasets/data/interim-full-moda>`_
- `ERA5 <http://apps.ecmwf.int/data-catalogues/era5/?class=ea>`_
- `MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_

Command Line Options
####################

- ``-D X``, ``--directory X``: working data directory
- ``-Y X``, ``--year X``: Years of model outputs to run
- ``--mean X``: Start and end year for mean
- ``--redistribute``: Uniformly redistribute values over the ocean
- ``-l X``, ``--lmax X``: maximum spherical harmonic degree
- ``-m X``, ``--mmax X``: maximum spherical harmonic order
- ``-n X``, ``--love X``: Load Love numbers dataset

    * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
    * ``1``: Gegout (2005) values from PREM [Gegout2010]_
    * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``-r X``, ``--reference X``: Reference frame for load love numbers

    * ``'CF'``: Center of Surface Figure (default)
    * ``'CM'``: Center of Mass of Earth System
    * ``'CE'``: Center of Mass of Solid Earth
- ``-F X``, ``--format X``: output data format

    * ``'ascii'``
    * ``'netCDF4'``
    * ``'HDF5'``
- ``-V``, ``--verbose``:  Output information for each output file
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Boy2005] J.-P. Boy and B. F. Chao, "Precise evaluation of atmospheric loading effects on Earth's time‐variable gravity field", *Journal of Geophysical Research: Solid Earth*, 110(B08412), (2005). `doi: 10.1029/2002JB002333 <https://doi.org/10.1029/2002JB002333>`_

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Swenson2002] S. Swenson and J. Wahr, "Estimated effects of the vertical structure of atmospheric mass on the time‐variable geoid", *Journal of Geophysical Research*, 107(B9), 2194, (2002) `doi: 10.1029/2000JB000024 <https://doi.org/10.1029/2000JB000024>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
