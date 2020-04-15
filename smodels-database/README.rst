================
SModelS Database
================

`SModelS`_ stores all the information about the experimental results in the Database.
The Database is organized as files in an ordinary (UNIX) directory hierarchy, with a thin Python
layer serving as the access to the database.
The overall structure of the directory hierarchy and its contents is
depicted in the scheme below:

.. image:: https://github.com/SModelS/smodels/blob/master/docs/manual/source/images/DatabaseFolders.png
   :scale: 30


The complete list of analyses and results included in the database can be
consulted at `https://smodels.github.io/docs/ListOfAnalyses <https://smodels.github.io/docs/ListOfAnalyses>`_.
We note that all the results in the official database release have been
carefully validated  and the validation material can be
found at `https://smodels.github.io/docs/Validation <https://smodels.github.io/docs/Validation>`_.


Installation
============

The database can conveniently be updated independently from `SModelS`_ code
updates. It suffices to download or clone this repository to a local folder and
correctly set the SModelS database path when running SModelS.
Alternatively, from `SModelS v1.1.3 <https://github.com/SModelS/smodels/releases>`_ onwards, the database path
can be specified as an URL, e.g. https://smodels.github.io/database/official123, and the binary
database file will be automatically downloaded and used. This is often faster than
building the binary file from the database folder and avoids possible machine dependences.
The database URLs can be found in the `releases page <https://github.com/SModelS/smodels-database-release/releases>`_.
For more information check the `SModelS online manual`_.


Adding FastLim data
^^^^^^^^^^^^^^^^^^^

The official SModelS database can be augmented with data from the
`fastlim <http://cern.ch/fastlim>`_ results.
For using SModelS with the text database,
a tarball with the *properly converted* fastlim-1.0 efficiency maps can be found in
the smodels-database folder.
The tarball then needs to be exploded in the top level directory of the database: ::

 cd <smodels-database folder>
 tar -xzvf smodels-v1.1-fastlim-1.0.tgz
 rm smodels-v1.1-fastlim-1.0.tgz

Once the fastlim folders have been added to the database,
SModelS auto-detects fastlim results and issues an acknowledgement.

As mentioned above, from `SModelS v1.1.3 <https://github.com/SModelS/smodels/releases>`_ onwards it is also possible to
directly download the database binary file using the URLs
provided in the `releases page <https://github.com/SModelS/smodels-database-release/releases>`_ .
Separate URLs are provided for the database including the Fastlim maps, so the user
can choose which database to use.

When using fastlim results, please properly cite the fastlim paper; for
convenience, a bibtex file is provided in the smodels-fastlim tarball.


Adding one's own results
^^^^^^^^^^^^^^^^^^^^^^^^

Adding additional experimental results is a matter of copying and editing text
files. Once the new folders and files have been added following the
database structure format, SModelS
automatically rebuilds the binary (Pickle) database file.
The added results will then be available for use with SModelS.


For citing the experimental analyses in the database, you can use
*database.bib*.

.. _SModelS online manual: https://smodels.readthedocs.io/
.. _SModelS: https://github.com/SModelS/smodels

