#!/usr/bin/env python3

"""
.. module:: testServer
   :synopsis: Tests the database server

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
import glob
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from smodels.experiment.databaseObj import Database
from unitTestHelpers import equalObjs, runMain, importModule
import time
import subprocess

from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.tools.databaseClient import DatabaseClient
from smodels.tools.databaseServer import DatabaseServer
        
setLogLevel ( "info" )

class ServerTest(unittest.TestCase):
    def testRun(self):
        port = 31770
        #server = DatabaseServer ( dbpath = "database/db31.pcl", port = port )
        #server.run ( nonblocking = True )
        startserver = f"../smodels/tools/smodelsTools.py proxydb -p {port} -i database/db31.pcl -o ./proxy31.pcl -r"
        cmd = startserver.split(" ")
        print ( "starting server %s" % startserver )
        subprocess.Popen ( cmd )
        print ( "started server" )
        time.sleep(3)
        filename = "./testFiles/slha/simplyGluino.slha"
        db = Database("./proxy31.pcl" )
        outputfile = runMain(filename,suppressStdout = False, overridedatabase = db )
        client = DatabaseClient ( port = port ) 
        client.send_shutdown()

if __name__ == "__main__":
    unittest.main()
