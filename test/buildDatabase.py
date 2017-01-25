#!/usr/bin/python3

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database

# dir = "./database/"
dir = "../../smodels-database/"

pcl = "%sdatabase.pcl" % dir

if os.path.exists ( pcl ):
    os.unlink ( pcl )

t0=time.time()
d=Database( dir )
print(d)
t1=time.time()
print ( "Building the database took %.2f seconds." % ( t1 - t0 ) )
s = os.stat ( pcl )
print ( "Database is %.1f MB." % ( s.st_size / 1000. / 1000. ) )
d=Database( dir )
t2=time.time()
print ( "Reading the database took %.2f seconds." % ( t2 - t1 ) )
