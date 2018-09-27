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
ver = "".join ( map ( str, version(True)[:3] ) )
dbname="http://smodels.hephy.at/database/unittest%s" % ver
database = Database(dbname, discard_zeroes = False)

if __name__ == "__main__":
    print ( database )
