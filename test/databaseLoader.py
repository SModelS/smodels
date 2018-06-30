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
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel('error')
database = Database( "./database", discard_zeroes = False, progressbar=True, force_load='txt')
