model-harmonics
===============

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/model-harmonics/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/model-harmonics/badge/?version=latest)](https://read-grace-harmonics.readthedocs.io/projects/model-harmonics/en/latest/?badge=latest)

Python tools for obtaining and working with model synthetic spherical harmonic coefficients for comparing with data from the the NASA/DLR Gravity Recovery and Climate Experiment (GRACE) and the NASA/GFZ Gravity Recovery and Climate Experiment Follow-On (GRACE-FO) missions

These are extension routines for the set of [read-GRACE-harmonics](https://github.com/tsutterley/read-GRACE-harmonics) tools

#### Resources
- [NASA GRACE mission site](https://www.nasa.gov/mission_pages/Grace/index.html)
- [NASA GRACE-FO mission site](https://www.nasa.gov/missions/grace-fo)
- [JPL GRACE Tellus site](https://grace.jpl.nasa.gov/)
- [JPL GRACE-FO site](https://gracefo.jpl.nasa.gov/)
- [UTCSR GRACE site](http://www.csr.utexas.edu/grace/)
- [GRACE at the NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)](https://podaac.jpl.nasa.gov/grace)
- [GRACE at the GFZ Information System and Data Center](http://isdc.gfz-potsdam.de/grace-isdc/)

#### Dependencies
- [geoid-toolkit: Python utilities for calculating geoid heights from static gravity field coefficients](https://github.com/tsutterley/geoid-toolkit/)
- [read-GRACE-harmonics: Python tools for working with GRACE/GRACE-FO data](https://github.com/tsutterley/read-GRACE-harmonics/)
- [pygrib: Python interface for reading and writing GRIB data](https://pypi.python.org/pypi/pygrib)
- [cdsapi: Python client libraries for the CDS Web API](https://pypi.org/project/cdsapi/)
- [fiona: Python wrapper for vector data access functions from the OGR library](https://fiona.readthedocs.io/en/latest/manual.html)
- [pyproj: Python interface to PROJ library](https://pypi.org/project/pyproj/)
- [shapely: PostGIS-ish operations outside a database context for Python](http://toblerity.org/shapely/index.html)
- [ecmwf-api-client: Python client libraries for the ECMWF Web API](https://software.ecmwf.int/wiki/display/WEBAPI/Web-API+Downloads)

#### References
I. Velicogna, Y. Mohajerani, G. A, F. Landerer, J. Mouginot, B. No&euml;l,
E. Rignot, T. C. Sutterley, M. van den Broeke, J. M. van Wessem, and D. Wiese,
"Continuity of ice sheet mass loss in Greenland and Antarctica from the GRACE
and GRACE Follow‐On missions", *Geophysical Research Letters*, 47,
(2020). [doi:10.1029/2020GL087291]( https://doi.org/10.1029/2020GL087291)

T. C. Sutterley, I. Velicogna, and C.-W. Hsu, "Self‐Consistent Ice Mass Balance
and Regional Sea Level From Time‐Variable Gravity", *Earth and Space Science*, 7,
(2020). [doi:10.1029/2019EA000860](https://doi.org/10.1029/2019EA000860)

#### Download
The program homepage is:  
https://github.com/tsutterley/model-harmonics  
A zip archive of the latest version is available directly at:  
https://github.com/tsutterley/model-harmonics/archive/main.zip  

#### Disclaimer
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the Center for Space Research at the University of Texas (UTCSR), the Jet Propulsion Laboratory (JPL), the German Research Centre for Geosciences (GeoForschungsZentrum, GFZ) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.

#### License
The content of this project is licensed under the [Creative Commons Attribution 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/) and the source code is licensed under the [MIT license](LICENSE).
