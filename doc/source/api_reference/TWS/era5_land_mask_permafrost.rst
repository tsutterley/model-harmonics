============================
era5_land_mask_permafrost.py
============================

- Creates a mask for ERA5-Land data based on the permafrost/surface classification from the `NSIDC Circum-Arctic Map of Permafrost and Ground-Ice Conditions <http://nsidc.org/data/ggd318.html>`_

    1. Continuous Permafrost
    2. Discontinuous Permafrost
    3. Isolated Permafrost
    4. Sporadic Permafrost
    5. Glaciated Area

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/TWS/era5_land_mask_permafrost.py

Calling Sequence
################

.. argparse::
    :filename: era5_land_mask_permafrost.py
    :func: arguments
    :prog: era5_land_mask_permafrost.py
    :nodescription:
    :nodefault:
