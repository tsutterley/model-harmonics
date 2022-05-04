=====================
ecco_read_realtime.py
=====================

- Reads 12-hour `ECCO ocean bottom pressure data from JPL <https://ecco.jpl.nasa.gov/drive/files/NearRealTime/Readme>`_
- Global area average of each ocean bottom pressure map is removed [Greatbatch1994]_
- Calculates monthly anomalies on an equirectangular grid

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/ECCO/ecco_read_realtime.py

Calling Sequence
################

.. argparse::
    :filename: ../ECCO/ecco_read_realtime.py
    :func: arguments
    :prog: ecco_read_realtime.py
    :nodescription:
    :nodefault:

    model : @after
        * ``'kf080i'``: Kalman filter analysis
        * ``'dr080i'``: RTS smoother analysis

References
##########

.. [Greatbatch1994] R. J. Greatbatch, "A note on the representation of steric sea level in models that conserve volume rather than mass", *Journal of Geophysical Research*, 99(C6), 12767--12771, (1994). `doi: 10.1029/94JC00847 <https://doi.org/10.1029/94JC00847>`_
