================
SModelS Database
================

This repository stores the private and development versions of the SModelS database.

Instructions for Releases of New Database Versions:
===================================================

#. Create a clean (public) version of the database
#. Merge or overwrite the latest version to the master branch of the `smodels-database-release repository <https://github.com/SModelS/smodels-database-release>`_
#. Go to the `smodels:master <https://github.com/SModelS/smodels/tree/master>`_ branch and run::

.. code-block::

   git subtree pull --prefix=smodels-database --squash git@github.com:SModelS/smodels-database-release.git master
   git push
   
The above will update the smodels-database folder (subtree) in smodels:master.   
