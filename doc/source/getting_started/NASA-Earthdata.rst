==============
NASA Earthdata
==============

NASA Data Distribution Centers
##############################

The NASA Earth Science Data Information Systems Project funds and operates
`12 Distributed Active Archive Centers (DAACs) <https://earthdata.nasa.gov/about/daacs>`_ throughout the United States.
These centers have recently transitioned from ftp to https servers.
The https updates are designed to increase performance and improve security during data retrieval.
NASA Earthdata uses `OAuth2 <https://wiki.earthdata.nasa.gov/pages/viewpage.action?pageId=71700485>`_, an approach to authentication that protects your personal information.

- https://urs.earthdata.nasa.gov/documentation
- https://wiki.earthdata.nasa.gov/display/EL/Knowledge+Base

GES DISC
########
The `Goddard Earth Sciences Data and Information Server Center (GES DISC) <https://disc.gsfc.nasa.gov/>`_ provides access to a wide range of global climate data, concentrated primarily in the areas of atmospheric composition, atmospheric dynamics, global precipitation, and solar irradiance.
If any problems contact GES DISC Help Desk at `gsfc-help-disc@lists.nasa.gov <mailto:gsfc-help-disc@lists.nasa.gov>`_ or the NASA EOSDIS support team `support@earthdata.nasa.gov <mailto:support@earthdata.nasa.gov>`_.
GES DISC support requests that you `include as much of the information as possible <https://disc.gsfc.nasa.gov/information/documents?title=Contact%20Us#email>`_ in your contact email.


ECCO Drive via PO.DAAC
######################

The `Physical Oceanography Distributed Active Archive Center (PO.DAAC) <https://podaac.jpl.nasa.gov/>`_ provides data and related information pertaining to the physical processes and conditions of the global oceans, including measurements of ocean winds, temperature, topography, salinity, circulation and currents, and sea ice.
If any problems contact JPL PO.DAAC support at `podaac@podaac.jpl.nasa.gov <mailto:podaac@podaac.jpl.nasa.gov>`_ or the NASA EOSDIS support team `support@earthdata.nasa.gov <mailto:support@earthdata.nasa.gov>`_.

WebDAV
------
JPL ECCO Drive uses passwords generated using the Web Distributed Authoring and Versioning (WebDAV) API.
This password is created at the `ECCO Drive <https://ecco.jpl.nasa.gov/drive>`_ website.
Use this password rather than your Earthdata password when retrieving model data from ECCO Drive.
`Click here for more information <https://ecco.jpl.nasa.gov/drive/help>`_.

Steps to Sync from GES DISC and PO.DAAC
#######################################

1. `Register with NASA Earthdata Login system <https://urs.earthdata.nasa.gov/users/new>`_
2. `After registering, login to the system <https://urs.earthdata.nasa.gov/home>`_
3. Add ``NASA GESDISC DATA ARCHIVE`` and ``ECCO Drive`` `applications to Earthdata <https://wiki.earthdata.nasa.gov/display/EL/How+To+Pre-authorize+an+application>`_
4. Sync `GLDAS land surface model outputs <https://github.com/tsutterley/model-harmonics/blob/main/GLDAS/gesdisc_gldas_sync.py>`_.
5. Retrieve `WebDAV password <https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_webdav.py>`_ to access ECCO servers
6. Sync `near real-time <https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_sync.py>`_ or `version 4 <https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_v4_sync.py>`_ data

Can also create a ``.netrc`` file for permanently storing NASA Earthdata and JPL WebDAV credentials:

.. code-block:: bash

    echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> ~/.netrc
    echo "machine ecco.jpl.nasa.gov login <uid> password <webdav>" >> ~/.netrc
    chmod 0600 ~/.netrc

Or set environmental variables for your NASA Earthdata and JPL WebDAV credentials:

.. code-block:: bash

    export EARTHDATA_USERNAME=<uid>
    export EARTHDATA_PASSWORD=<password>
    export ECCO_PASSWORD=<webdav>

Other Data Access Examples
##########################
- `Curl and Wget <https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget>`_
- `Python <https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python>`_
