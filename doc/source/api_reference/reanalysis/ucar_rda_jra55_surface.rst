=========================
ucar_rda_jra55_surface.py
=========================

- Downloads JRA-55 products using a links list csh file provided by the `NCAR/UCAR Research Data Archive (RDA) <https://rda.ucar.edu/>`_
- JRA-55 6-hour data is more regularly updated compared with the monthly means
- Combines 6-hour model outputs into monthly averages
- Will extract data files if compressed into a single tar file

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/reanalysis/ucar_rda_jra55_surface.py

Calling Sequence
################

.. argparse::
    :filename: ucar_rda_jra55_surface.py
    :func: arguments
    :prog: ucar_rda_jra55_surface.py
    :nodescription:
    :nodefault:

Notes
#####
- Need to make small enough requests so that the data comes in a single tar file or as single files
- Large requests coming in multiple tar files are presently not supported
