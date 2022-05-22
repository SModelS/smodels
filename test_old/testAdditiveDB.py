#!/usr/bin/env python3

"""
.. module:: testAdditiveDB
   :synopsis: performs tests with adding databases

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
import unittest
import logging.config
import os

class AdditiveDBTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    from smodels.tools.smodelsLogging import logger

    def testAddingSubset(self):
        """ tests adding same dataset """
        db1 = Database ( "./database/" )
        db2 = Database ( "./tinydb/+./database/" ) ## tinydb is subset of database
        self.assertTrue ( len(db2.expResultList) == len(db1.expResultList) )
        tx1, tx2 = [], []
        for er1, er2 in zip ( db1.getExpResults(), db2.getExpResults() ):
            tx1.append ( er1.getTxNames() )
            tx2.append ( er2.getTxNames() )
        self.assertTrue ( len(tx1) == len(tx2) )

    def testAddingTxname(self):
        """ tests adding same dataset """
        sel = ["CMS-PAS-SUS-15-002" ]
        tx1, tx2 = [], []
        db1 = Database ( "./database/" )
        for er1 in db1.getExpResults( analysisIDs= sel ): 
            tx1.append ( er1.getTxNames() )
        tx1 = tx1[0]
        db2 = Database ( "./database/+./dbadd1/" ) ## adds a t1 txname
        self.assertTrue ( len(db2.expResultList) == len(db1.expResultList) )
        for er2 in db2.getExpResults( analysisIDs = sel ):
            tx2.append ( er2.getTxNames() )
        tx2 = tx2[0]
        self.assertTrue ( len(tx2) == len(tx1)+1 )

    def testAddingDataset(self):
        """ tests adding a dataset """
        sel = ["CMS-SUS-13-012" ]
        tx1, tx2 = [], []
        db1 = Database ( "./database/" )
        for er1 in db1.getExpResults( analysisIDs= sel ): 
            tx1.append ( er1.getTxNames() )
        tx1 = tx1[0]
        db2 = Database ( "./database/+./dbadd1/" ) ## adds a t1 txname
        self.assertTrue ( len(db2.expResultList) == len(db1.expResultList) )
        for er2 in db2.getExpResults( analysisIDs = sel ):
            tx2.append ( er2.getTxNames() )
        tx2 = tx2[0]
        db2 = Database ( "./database/+./dbadd1/" ) ## adds a t1 txname
        self.assertTrue ( len(tx2) == len(tx1)+5 )

if __name__ == "__main__":
    unittest.main()
