======================
ecco_read_llc_tiles.py
======================

- Reads monthly ECCO ocean bottom pressure LLC tile data and calculates multi-annual means
- Global area average of each ocean bottom pressure map is removed [Greatbatch1994]_

Calling Sequence

.. code-block:: bash

    python ecco_read_llc_tiles.py --directory <path_to_directory> V4r4 V5alpha

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_llc_tiles.py


Inputs
######

- ECCO Version 4 or 5 models

    * ``'V4r4'``: Version 4, Revision 4
    * ``'V5alpha'``: ECCO Version 5, Alpha release

Command Line Options
####################

- ``-D X``, ``--directory X``: working data directory
- ``-Y X``, ``--year X``: Years to run
- ``-m X``, ``--mean X``: Year range for mean
- ``-M X``, ``--mode X``: Permission mode of directories and files
- ``-V``, ``--verbose``: Output information for each output file

References
##########

.. [Greatbatch1994] R. J. Greatbatch, "A note on the representation of steric sea level in models that conserve volume rather than mass", *Journal of Geophysical Research*, 99(C6), 12767--12771, (1994). `doi: 10.1029/94JC00847 <https://doi.org/10.1029/94JC00847>`_
