===================
GIA_GSFC_mascons.py
===================

- Calculates GIA equivalent water height corrections for GSFC mascons at the central points of each mascon

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/GIA/GIA_GSFC_mascons.py

Calling Sequence
################

.. argparse::
    :filename: GIA_GSFC_mascons.py
    :func: arguments
    :prog: GIA_GSFC_mascons.py
    :nodescription:
    :nodefault:

    --gia -G : @after
        * ``'IJ05-R2'``: `Ivins R2 GIA Models <https://doi.org/10.1002/jgrb.50208>`_
        * ``'W12a'``: `Whitehouse GIA Models <https://doi.org/10.1111/j.1365-246X.2012.05557.x>`_
        * ``'SM09'``: `Simpson/Milne GIA Models <https://doi.org/10.1029/2010JB007776>`_
        * ``'ICE6G'``: `ICE-6G GIA Models <https://doi.org/10.1002/2014JB011176>`_
        * ``'Wu10'``: `Wu (2010) GIA Correction <https://doi.org/10.1038/ngeo938>`_
        * ``'AW13-ICE6G'``: `Geruo A ICE-6G GIA Models <https://doi.org/10.1093/gji/ggs030>`_
        * ``'AW13-IJ05'``: `Geruo A IJ05-R2 GIA Models <https://doi.org/10.1093/gji/ggs030>`_
        * ``'Caron'``: `Caron JPL GIA Assimilation <https://doi.org/10.1002/2017GL076644>`_
        * ``'ICE6G-D'``: `ICE-6G Version-D GIA Models <https://doi.org/10.1002/2016JB013844>`_
        * ``'ascii'``: reformatted GIA in ascii format
        * ``'netCDF4'``: reformatted GIA in netCDF4 format
        * ``'HDF5'``: reformatted GIA in HDF5 format
