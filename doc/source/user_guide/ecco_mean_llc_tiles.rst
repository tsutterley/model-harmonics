======================
ecco_mean_llc_tiles.py
======================

- Reads monthly ECCO ocean bottom pressure LLC tile data and calculates multi-annual means
- Global area average of each ocean bottom pressure map is removed [Greatbatch1994]_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_mean_llc_tiles.py

Calling Sequence
################

.. argparse::
    :filename: ../ECCO/ecco_mean_llc_tiles.py
    :func: arguments
    :prog: ecco_mean_llc_tiles.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'V4r4'``: Version 4, Revision 4
        * ``'V5alpha'``: ECCO Version 5, Alpha release

References
##########

.. [Greatbatch1994] R. J. Greatbatch, "A note on the representation of steric sea level in models that conserve volume rather than mass", *Journal of Geophysical Research*, 99(C6), 12767--12771, (1994). `doi: 10.1029/94JC00847 <https://doi.org/10.1029/94JC00847>`_
