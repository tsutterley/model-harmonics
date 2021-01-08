============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Adds additional modules to the GRACE/GRACE-FO `gravity_toolkit utilities <https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/utilities.py>`__

 - Can list a directory from the `Goddard Earth Sciences Data and Information Server Center (GES DISC) <https://disc.gsfc.nasa.gov/>`_

`Source code`__

.. __: https://github.com/tsutterley/model-harmonics/blob/main/model_toolkit/utilities.py


General Methods
===============


.. method:: gravity_toolkit.utilities.gesdisc_list(HOST,username=None,password=None,build=True,timeout=None,urs=None,parser=None,format=None,pattern='',sort=False)

    List a directory on NASA `GES DISC <https://disc.gsfc.nasa.gov/>`_ servers

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `build`: Build opener with NASA Earthdata credentials

        `urs`: Earthdata login URS 3 host

        `timeout`: timeout in seconds for blocking operations

        `parser`: HTML parser for lxml

        `format`: format for input time string

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `colnames`: list of column names in a directory

        `collastmod`: list of last modification times for items in the directory



