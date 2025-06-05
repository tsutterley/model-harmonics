===============
model-harmonics
===============

|License|
|Documentation Status|
|PyPI|
|commits-since|
|zenodo|

.. |License| image:: https://img.shields.io/github/license/tsutterley/model-harmonics
   :target: https://github.com/tsutterley/model-harmonics/blob/main/LICENSE

.. |Documentation Status| image:: https://readthedocs.org/projects/model-harmonics/badge/?version=latest
   :target: https://model-harmonics.readthedocs.io/en/latest/?badge=latest

.. |PyPI| image:: https://img.shields.io/pypi/v/model-harmonics.svg
   :target: https://pypi.python.org/pypi/model-harmonics/

.. |commits-since| image:: https://img.shields.io/github/commits-since/tsutterley/model-harmonics/latest
   :target: https://github.com/tsutterley/model-harmonics/releases/latest

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5156964.svg
   :target: https://doi.org/10.5281/zenodo.5156964


Python tools for obtaining and working with model synthetic spherical harmonic coefficients for comparing with data from the the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

These are extension routines for the set of `gravity-toolkit <https://github.com/tsutterley/gravity-toolkit>`_ tools

Resources
#########

- `NASA GRACE mission site <https://www.nasa.gov/mission_pages/Grace/index.html>`_
- `NASA GRACE-FO mission site <https://www.nasa.gov/missions/grace-fo>`_
- `JPL GRACE Tellus site <https://grace.jpl.nasa.gov/>`_
- `JPL GRACE-FO site <https://gracefo.jpl.nasa.gov/>`_
- `UTCSR GRACE site <http://www.csr.utexas.edu/grace/>`_
- `GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC) <https://podaac.jpl.nasa.gov/grace>`_
- `GRACE at the GFZ Information System and Data Center <http://isdc.gfz-potsdam.de/grace-isdc/>`_

Dependencies
############

- `geoid-toolkit: Python utilities for calculating geoid heights from static gravity field coefficients <https://github.com/tsutterley/geoid-toolkit/>`_
- `gravity-toolkit: Python tools for working with GRACE/GRACE-FO data <https://github.com/tsutterley/gravity-toolkit/>`_
- `pygrib: Python interface for reading and writing GRIB data <https://pypi.python.org/pypi/pygrib>`_
- `cdsapi: Python client libraries for the CDS Web API <https://pypi.org/project/cdsapi/>`_
- `fiona: Python wrapper for vector data access functions from the OGR library <https://fiona.readthedocs.io/en/latest/manual.html>`_
- `pyproj: Python interface to PROJ library <https://pypi.org/project/pyproj/>`_
- `shapely: PostGIS-ish operations outside a database context for Python <http://toblerity.org/shapely/index.html>`_
- `ecmwf-api-client: Python client libraries for the ECMWF Web API <https://software.ecmwf.int/wiki/display/WEBAPI/Web-API+Downloads>`_

References
##########

    I. Velicogna, Y. Mohajerani, G. A, F. Landerer, J. Mouginot, B. No&euml;l,
    E. Rignot, T. C. Sutterley, M. van den Broeke, J. M. van Wessem, and D. Wiese,
    "Continuity of ice sheet mass loss in Greenland and Antarctica from the GRACE
    and GRACE Follow‐On missions", *Geophysical Research Letters*, 47,
    (2020). `doi: 10.1029/2020GL087291 <https://doi.org/10.1029/2020GL087291>`_

    T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self‐Consistent Ice Mass Balance
    and Regional Sea Level From Time‐Variable Gravity", *Earth and Space Science*, 7,
    (2020). `doi: 10.1029/2019EA000860 <https://doi.org/10.1029/2019EA000860>`_

Download
########

| The program homepage is:
| https://github.com/tsutterley/model-harmonics
| A zip archive of the latest version is available directly at:
| https://github.com/tsutterley/model-harmonics/archive/main.zip

Disclaimer
##########

This project contains work and contributions from the `scientific community <./CONTRIBUTORS.rst>`_.
This program is not sponsored or maintained by the Universities Space Research Association (USRA),
the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL),
the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.
It is provided here for your convenience but *with no guarantees whatsoever*.

License
#######

The content of this project is licensed under the `Creative Commons Attribution 4.0 Attribution license <https://creativecommons.org/licenses/by/4.0/>`_ and the source code is licensed under the `MIT license <LICENSE>`_.
