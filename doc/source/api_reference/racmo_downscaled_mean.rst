========================
racmo_downscaled_mean.py
========================

- Calculates the temporal mean of downscaled RACMO surface mass balance products

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/SMB/racmo_downscaled_mean.py

Calling Sequence
################

.. argparse::
    :filename: racmo_downscaled_mean.py
    :func: arguments
    :prog: racmo_downscaled_mean.py
    :nodescription:
    :nodefault:

    --product -P : @after
        * ``'PRECIP'``: Precipitation
        * ``'REFREEZE'``: Meltwater Refreeze
        * ``'RUNOFF'``: Meltwater Runoff
        * ``'SMB'``: Surface Mass Balance
        * ``'SNOWMELT'``: Snowmelt
