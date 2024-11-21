===============
regress_maps.py
===============

- Reads in spatial files and fits a regression model at each grid point

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/scripts/regress_maps.py

Calling Sequence
################

.. argparse::
    :filename: regress_maps.py
    :func: arguments
    :prog: regress_maps.py
    :nodescription:
    :nodefault:

    --interval : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)
