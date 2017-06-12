#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
colors.on = True
setLogLevel ( "debug" )

# dir = "../../smodels-database/"
dir = "tinydb/"
d=Database( dir, discard_zeroes = True )

# print ( "version", d.databaseVersion )
