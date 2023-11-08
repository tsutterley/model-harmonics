================
jpl_ecco_sync.py
================

- Syncs ECCO Near Real-Time model outputs from the `NASA JPL ECCO Drive server <https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme>`_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_sync.py

Calling Sequence
################

.. argparse::
    :filename: jpl_ecco_sync.py
    :func: arguments
    :prog: jpl_ecco_sync.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'kf080i'``: Kalman filter analysis
        * ``'dr080i'``: RTS smoother analysis
