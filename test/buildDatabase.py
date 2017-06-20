#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
colors.on = True
setLogLevel ( "debug" )

dir = "corrdb/"
d=Database( dir, discard_zeroes = True )
print(d)
e=d.getExpResults()
e0=e[0]
print (e0)
print ( "upper limit", e0.getUpperLimitFor ( dataID=None ) )
