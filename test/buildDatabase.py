#!/usr/bin/python

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel ( "debug" )

# d=Database("./database/database%d.pcl" % int ( sys.version[0] ) )
d=Database("./database/" )
print(d)
