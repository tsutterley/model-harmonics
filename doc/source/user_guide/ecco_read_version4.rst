=====================
ecco_read_version4.py
=====================

- Reads monthly ECCO ocean bottom pressure data from `Version 4 models <https://ecco-group.org/products-ECCO-V4r4.htm>`_ and calculates monthly anomalies
- Global area average of each ocean bottom pressure map is removed [Greatbatch1994]_

Calling Sequence

.. code-block:: bash

    python ecco_read_version4.py --directory <path_to_directory> V4r3 V4r4

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_version4.py

Inputs
######

- ECCO Version 4 models

    * ``'V4r3'``: Version 4, Revision 3
    * ``'V4r4'``: Version 4, Revision 4

Command Line Options
####################

- ``-D X``, ``--directory X``: working data directory
- ``-Y X``, ``--year X``: Years to run
- ``-m X``, ``--mean X``: Year range for mean
- ``-F X``, ``--format X``: input and output data format

    * ``'ascii'``
    * ``'netcdf'``
    * ``'HDF5'``
- ``-M X``, ``--mode X``: Permission mode of directories and files
- ``-V``, ``--verbose``: Output information for each output file

References
##########

.. [Greatbatch1994] R. J. Greatbatch, "A note on the representation of steric sea level in models that conserve volume rather than mass", *Journal of Geophysical Research*, 99(C6), 12767--12771, (1994). `doi: 10.1029/94JC00847 <https://doi.org/10.1029/94JC00847>`_
