==================================
``reanalysis_monthly_pressure.py``
==================================

Reads daily atmospheric pressure fields from reanalysis and outputs monthly averages

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/model_harmonics/reanalysis/reanalysis_monthly_pressure.py

Calling Sequence
################

.. argparse::
    :module: model_harmonics.reanalysis.reanalysis_mean_pressure
    :func: arguments
    :prog: reanalysis_mean_pressure.py
    :nodescription:
    :nodefault:

    model : @after
        * `NCEP-DOE-2 <https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.html>`_
