#!/usr/bin/env python

"""
.. module:: databaseLoader
   :synopsis: When running the complete test suite, we need to
              load the database only once

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
database = Database ( "./database/db%d0.pcl" % int ( sys.version[0] ), \
                      discard_zeroes = False  )
