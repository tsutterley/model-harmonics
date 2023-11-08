==============================
racmo_downscaled_cumulative.py
==============================

- Calculates cumulative anomalies of downscaled RACMO surface mass balance products

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/SMB/racmo_downscaled_cumulative.py

Calling Sequence
################

.. argparse::
    :filename: racmo_downscaled_cumulative.py
    :func: arguments
    :prog: racmo_downscaled_cumulative.py
    :nodescription:
    :nodefault:

    --product -P : @after
        * ``'PRECIP'``: Precipitation
        * ``'REFREEZE'``: Meltwater Refreeze
        * ``'RUNOFF'``: Meltwater Runoff
        * ``'SMB'``: Surface Mass Balance
        * ``'SNOWMELT'``: Snowmelt
