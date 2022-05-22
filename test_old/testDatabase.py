#!/usr/bin/env python3

"""
.. module:: testDatabase
   :synopsis: performs tests with database loading, pickle writing, filtering

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
import unittest
import logging.config
import os


class DatabaseTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    from smodels.tools.smodelsLogging import logger

    def testWritePickle(self):
        """ tests writing pickle file """
        binfile = "./.database.pcl"
        if os.path.exists ( binfile ):
            os.unlink ( binfile )
        self.logger.info ( "test writing pickle file """ )
        writer = Database ( "./tinydb/", force_load = "txt" )
        writer.createBinaryFile ( binfile )
        reader = Database ( binfile, force_load="pcl" )

        os.unlink ( binfile )
        self.assertEqual( writer, reader )

    def testSelectors(self):
        from databaseLoader import database
        validated = database.getExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'],
                                            useNonValidated = False )
        nonval = database.getExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'],
                                         useNonValidated = True )
#         print ( "validated=",len(validated),map ( str, validated ) )
        self.assertTrue(len(validated)==8)
        self.assertTrue(len(nonval)==9)


    def testLoadLatest(self):
        dblatest=Database ("latest")
        latestver = dblatest.databaseVersion.replace(".","").replace("official","")
        from databaseLoader import database
        thisver = database.databaseVersion.replace("unittest","").replace(".","")
        ilatestver = int ( latestver[:2] )
        ithisver = int ( thisver[:2] )
        if ilatestver < ithisver-1:
            print ( "latest and this version are different. is that an issue?",
                    ilatestver,"versus",ithisver )
        self.assertTrue ( ilatestver>=(ithisver-1) )


if __name__ == "__main__":
    unittest.main()
