=====================
ecco_read_version4.py
=====================

- Reads monthly ECCO ocean bottom pressure data from `Version 4 models <https://ecco-group.org/products-ECCO-V4r4.htm>`_ and calculates monthly anomalies
- Global area average of each ocean bottom pressure map is removed [Greatbatch1994]_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_version4.py

Calling Sequence
################

.. argparse::
    :filename: ../ECCO/ecco_read_version4.py
    :func: arguments
    :prog: ecco_read_version4.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'V4r3'``: Version 4, Revision 3
        * ``'V4r4'``: Version 4, Revision 4

References
##########

.. [Greatbatch1994] R. J. Greatbatch, "A note on the representation of steric sea level in models that conserve volume rather than mass", *Journal of Geophysical Research*, 99(C6), 12767--12771, (1994). `doi: 10.1029/94JC00847 <https://doi.org/10.1029/94JC00847>`_
