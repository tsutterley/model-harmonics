======================
merra_hybrid_regrid.py
======================

- Read and regrid MERRA-2 hybrid derived surface mass balance products

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_hybrid_regrid.py

Calling Sequence
################

.. argparse::
    :filename: merra_hybrid_regrid.py
    :func: arguments
    :prog: merra_hybrid_regrid.py
    :nodescription:
    :nodefault:

    --interval -I : @replace
        Output grid interval

        * ``1``: (0:360, 90:-90)
        * ``2``: (degree spacing/2)
        * ``3``: non-global grid (set with defined bounds)
