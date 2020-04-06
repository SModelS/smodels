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
# database = Database("./database", discard_zeroes = False)
# database = Database("../../smodels-database", discard_zeroes = False)
database = Database( "unittest", discard_zeroes = False)

if __name__ == "__main__":
    print ( database )
