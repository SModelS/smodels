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

# dbpath = "./database"
# dbpath = "../../smodels-database"
dbpath = "unittest"

database = Database( dbpath, discard_zeroes = False)

if __name__ == "__main__":
    print ( database )
