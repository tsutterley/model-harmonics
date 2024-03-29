========================
gldas_scaling_factors.py
========================

- Reads monthly GLDAS total water storage anomalies and monthly spherical harmonic coefficients [Wahr1998]_
- Computes point scaling factors following [Landerer2012]_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gldas_scaling_factors.py

Calling Sequence
################

.. argparse::
    :filename: gldas_scaling_factors.py
    :func: arguments
    :prog: gldas_scaling_factors.py
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

    --love -n : @after
        * ``0``: Han and Wahr (1995) values from PREM [Han1995]_
        * ``1``: Gegout (2005) values from PREM [Gegout2010]_
        * ``2``: Wang et al. (2012) values from PREM [Wang2012]_
        * ``3``: Wang et al. (2012) values from PREM with hard sediment [Wang2012]_
        * ``4``: Wang et al. (2012) values from PREM with soft sediment [Wang2012]_

    --reference : @after
        * ``'CF'``: Center of Surface Figure
        * ``'CM'``: Center of Mass of Earth System
        * ``'CE'``: Center of Mass of Solid Earth

References
##########

.. [Gegout2010] P. Gegout, J. Boehm, and D. Wijaya, "Practical numerical computation of love numbers and applications", Workshop of the COST Action ES0701, (2010). `doi: 10.13140/RG.2.1.1866.7045 <https://doi.org/10.13140/RG.2.1.1866.7045>`_

.. [Han1995] D. Han and J. Wahr, "The viscoelastic relaxation of a realistically stratified earth, and a further analysis of postglacial rebound", *Geophysical Journal International*, 120(2), 287--311, (1995). `doi: 10.1111/j.1365-246X.1995.tb01819.x <https://doi.org/10.1111/j.1365-246X.1995.tb01819.x>`_

.. [Landerer2012] F. W. Landerer and S. C. Swenson, "Accuracy of scaled GRACE terrestrial water storage estimates", *Water Resources Research*, 48(W04531), (2012). `doi: 10.1029/2011WR011453 <https://doi.org/10.1029/2011WR011453>`_

.. [Wahr1998] J. Wahr, M. Molenaar, and F. Bryan, "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE", *Journal of Geophysical Research*, 103(B12), 30205--30229, (1998). `doi: 10.1029/98JB02844 <https://doi.org/10.1029/98JB02844>`_

.. [Wang2012] H. Wang et al., "Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and modified models with refined crustal structure from Crust 2.0", *Computers & Geosciences*, 49, 190--199, (2012). `doi: 10.1016/j.cageo.2012.06.022 <https://doi.org/10.1016/j.cageo.2012.06.022>`_
