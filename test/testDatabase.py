#!/usr/bin/env python

"""
.. module:: testDatabase
   :synopsis: performs tests with database loading, pickle writing, ....
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
import unittest
import logging
import logging.config
import os

class DatabaseTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    logger = logging.getLogger(__name__)

    def testWritePickle(self):
        """ tests writing pickle file """
        binfile = "./.database.pcl"
        if os.path.exists ( binfile ):
            os.unlink ( binfile )
        self.logger.info ( "test writing pickle file """ )
        writer = Database ( "./database/", force_load = "txt" )
        writer.createBinaryFile ( binfile )
        reader = Database ( binfile )
        os.unlink ( binfile )
        self.assertEqual( writer, reader )

if __name__ == "__main__":
    unittest.main()
