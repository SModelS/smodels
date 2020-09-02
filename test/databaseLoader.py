#!/usr/bin/env python3

"""
.. module:: databaseLoader
   :synopsis: When running the complete test suite, we need to
              load the database only once

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.installation import version
ver = version().replace(".","")
#dbname="./database/db%d0.pcl" % int ( sys.version[0] )
dbname="https://smodels.github.io/database/unittest%s" % ver
#dbname="database/"
database = Database(dbname, discard_zeroes = False)

if __name__ == "__main__":
    print ( database )
