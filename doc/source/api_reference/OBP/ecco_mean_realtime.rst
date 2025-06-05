=====================
ecco_mean_realtime.py
=====================

- Reads 12-hour `ECCO ocean bottom pressure data from JPL <https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme>`_
- Global area average of each ocean bottom pressure map is removed :cite:p:`Greatbatch:1994dx`
- Calculates multi-annual means on an equirectangular grid

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/OBP/ecco_mean_realtime.py

Calling Sequence
################

.. argparse::
    :filename: ecco_mean_realtime.py
    :func: arguments
    :prog: ecco_mean_realtime.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'kf080i'``: Kalman filter analysis
        * ``'dr080i'``: RTS smoother analysis
