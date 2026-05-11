#!/usr/bin/env python3

"""
.. module:: fullDatabase
   :synopsis: Some unit tests require the full database, so lets give
              a centralized version

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database

dbpath = "official"
dbpath = "../../branches/smodels-database/"
database = Database( dbpath)

if __name__ == "__main__":
    print ( database )
