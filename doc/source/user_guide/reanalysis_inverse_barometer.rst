===============================
reanalysis_inverse_barometer.py
===============================

- Reads hourly mean sea level pressure fields from reanalysis and calculates the inverse-barometer response [Wunsch1997]_ [HofmannWellenhof2006]_

Calling Sequence
################

.. code-block:: bash

    python reanalysis_inverse_barometer.py --directory <path_to_directory> ERA5 MERRA-2

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_inverse_barometer.py

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
- ``-d X``, ``--density X``: Density of seawater in kg/m\ :sup:`3`
- ``-V``, ``--verbose``:  Output information for each output file
- ``-M X``, ``--mode X``: Permissions mode of the files created

References
##########

.. [Wunsch1997] Wunsch and Stammer. "Atmospheric loading and the oceanic "inverted barometer" effect", *Reviews of Geophysics*, 35(1), 79-107, (1997). `doi:10.1029/96RG03037 <https://doi.org/10.1029/96RG03037>`_
.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
