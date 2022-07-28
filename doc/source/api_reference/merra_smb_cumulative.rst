=======================
merra_smb_cumulative.py
=======================

- Reads MERRA-2 datafiles to calculate monthly cumulative anomalies in derived surface mass balance products

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/SMB/merra_smb_cumulative.py

Calling Sequence
################

.. argparse::
    :filename: ../../SMB/merra_smb_cumulative.py
    :func: arguments
    :prog: merra_smb_cumulative.py
    :nodescription:
    :nodefault:

    product : @after
        * ``'SMB'``: Surface Mass Balance
        * ``'ACCUM'``: Snowfall accumulation
        * ``'PRECIP'``: Total Precipitation
        * ``'RAINFALL'``: Total Rainfall
        * ``'SUBLIM'``: Evaporation and Sublimation
        * ``'RUNOFF'``: Meltwater Runoff
