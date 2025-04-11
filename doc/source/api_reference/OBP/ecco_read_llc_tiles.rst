======================
ecco_read_llc_tiles.py
======================

- Reads monthly ECCO ocean bottom pressure LLC tile data and calculates multi-annual means
- Global area average of each ocean bottom pressure map is removed :cite:p:`Greatbatch:1994dx`

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/OBP/ecco_read_llc_tiles.py

Calling Sequence
################

.. argparse::
    :filename: ecco_read_llc_tiles.py
    :func: arguments
    :prog: ecco_read_llc_tiles.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'V4r4'``: Version 4, Revision 4
        * ``'V5alpha'``: ECCO Version 5, Alpha release
