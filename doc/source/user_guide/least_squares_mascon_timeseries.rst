==================================
least_squares_mascon_timeseries.py
==================================

- Reads in an index of spherical harmonic coefficient files
- Filters and smooths data with specified processing algorithms [Jekeli1981]_ [Swenson2006]_
- Calculates a time-series of regional mass anomalies through a least-squares mascon procedure [Tiwari2009]_ [Jacob2012]_

Calling Sequence
################

.. code-block:: bash

    python least_squares_mascon_timeseries.py --mode 0o775 parameter_file

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/scripts/least_squares_mascon_timeseries.py

Inputs
######

- parameter file containing specific variables for the analysis

Command Line Options
####################

- ``-P X``, ``--np X``: Run in parallel with X number of processes
- ``-n X``, ``--love X``: Load Love numbers dataset

     * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
     * ``1``: Gegout (2005) values from PREM [Gegout2010]_
     * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
- ``-r X``, ``--reference X``: Reference frame for load love numbers

     * ``'CF'``: Center of Surface Figure (default)
     * ``'CM'``: Center of Mass of Earth System
     * ``'CE'``: Center of Mass of Solid Earth
- ``-l``, ``--log``: output log file for each job
- ``-M X``, ``--mode X``: permissions mode of output files

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Jacob2012] T. Jacob, J. Wahr, W. T. Pfeffer, and S. Swenson, "Recent contributions of glaciers and ice caps to sea level rise", *Nature*, 482, 514--518, (2012). `doi: 10.1038/nature10847 <https://doi.org/10.1038/nature10847>`_

.. [Jekeli1981] C. Jekeli, "Alternative Methods to Smooth the Earth's Gravity Field", NASA Grant No. NGR 36-008-161, OSURF Proj. No. 783210, 48 pp., (1981).

.. [Swenson2006] S. Swenson and J. Wahr, "Post‐processing removal of correlated errors in GRACE data", *Geophysical Research Letters*, 33(L08402), (2006). `doi: 10.1029/2005GL025285 <https://doi.org/10.1029/2005GL025285>`_

.. [Tiwari2009] V. M. Tiwari, J. Wahr, and S. Swenson, "Dwindling groundwater resources in northern India, from satellite gravity observations", *Geophysical Research Letters*, 36(L18401), (2009). `doi: 10.1029/2009GL039401 <https://doi.org/10.1029/2009GL039401>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
