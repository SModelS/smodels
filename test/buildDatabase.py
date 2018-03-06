#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
colors.on = True
setLogLevel ( "info" )

dir = "../../smodels-database/"
dir = "database/"
discard_zeroes = True

t0=time.time()
pcl = dir + "db%d%d.pcl" % ( int ( sys.version[0] ), int ( discard_zeroes ) )
if os.path.exists ( pcl ):
    os.unlink ( pcl )
d=Database( dir, discard_zeroes = discard_zeroes, progressbar = True )
print(d)
t1=time.time()
print ( "Building the database took %.2f seconds." % ( t1 - t0 ) )
s = os.stat ( pcl )
print ( "Database is %.1f MB." % ( s.st_size / 1024. / 1024. ) )
d=Database( pcl, discard_zeroes = True, force_load="pcl" )
t2=time.time()
print ( "Reading the database took %.2f seconds." % ( t2 - t1 ) )
