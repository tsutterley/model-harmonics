====================
spatial_operators.py
====================

- Performs basic operations on spatial files

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/scripts/spatial_operators.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/spatial_operators.py
    :func: arguments
    :prog: spatial_operators.py
    :nodescription:
    :nodefault:

    --interval -I : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
