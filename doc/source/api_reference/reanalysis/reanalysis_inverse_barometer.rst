===============================
reanalysis_inverse_barometer.py
===============================

- Reads hourly mean sea level pressure fields from reanalysis and calculates the inverse-barometer response :cite:p:`Wunsch:1997kg` :cite:p:`HofmannWellenhof:2006hy`

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/reanalysis_inverse_barometer.py

Calling Sequence
################

.. argparse::
    :filename: reanalysis_inverse_barometer.py
    :func: arguments
    :prog: reanalysis_inverse_barometer.py
    :nodescription:
    :nodefault:

    model : @after
        * `ERA-Interim <http://apps.ecmwf.int/datasets/data/interim-full-moda>`_
        * `ERA5 <http://apps.ecmwf.int/data-catalogues/era5/?class=ea>`_
        * `MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_

    --density -D : @replace
        Density of seawater in kg/m\ :sup:`3`
