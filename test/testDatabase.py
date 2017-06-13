#!/usr/bin/env python

"""
.. module:: testDatabase
   :synopsis: performs tests with database loading, pickle writing, filtering
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
import unittest
import logging.config
import os

class DatabaseTest(unittest.TestCase):
    # use different logging config for the unit tests.

    from smodels.tools.smodelsLogging import getLogger
    logger = getLogger ( toFile="smodels.log" )

    def testWritePickle(self):
        """ tests writing pickle file """
        binfile = "./.database.pcl"
        if os.path.exists ( binfile ):
            os.unlink ( binfile )
        self.logger.info ( "test writing pickle file """ )
        writer = Database ( "./tinydb/", force_load = "txt" )
        # writer = Database ( "./tinydb/" )
        writer.createBinaryFile ( binfile )
        reader = Database ( binfile, force_load="pcl" )
        os.unlink ( binfile )
        self.assertEqual( writer, reader )

    def testSelectors(self):
        from databaseLoader import database 
        validated = database.getExpResults ( useNonValidated = False )
        nonval = database.getExpResults ( useNonValidated = True )
        #print ( "validated=",len(validated),map ( str, validated ) )
        #print ( "non val=",len(nonval) )
        self.assertTrue ( len(validated)==7 )
        self.assertTrue ( len(nonval)==8 )





if __name__ == "__main__":
    unittest.main()
