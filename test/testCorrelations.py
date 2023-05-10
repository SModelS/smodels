#!/usr/bin/env python3

"""
.. module:: testCorrelations
   :synopsis: Test the correlation checks

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from databaseLoader import database

import unittest

class CorrelationsTest(unittest.TestCase):
    def testDifferentExperiments(self):
        ds1 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-05" ] )[0].datasets[0]
        ds2 = database.getExpResults( analysisIDs = [ "CMS-EXO-13-006" ] )[0].datasets[0]
        self.assertTrue ( ds1.isCombinableWith ( ds2 ) )

    def testDifferentRuns(self):
        ds1 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-05" ] )[0].datasets[0]
        ds2 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2016-08" ] )[0].datasets[0]
        self.assertTrue ( ds1.isCombinableWith ( ds2 ) )

    def testViaCombinationsMatrix(self):
        """ check that combinability via "combinableWith" fields works 
        in the test base, ATLAS-SUSY-2013-05 is marked as combinable with ATLAS-SUSY-2013-02
        """
        ds1 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-05" ] )[0].datasets[0]
        ds3 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-12" ] )[0].datasets[0]
        self.assertFalse ( ds1.isCombinableWith ( ds3 ) )
        database.combinationsmatrix = { "ATLAS-SUSY-2013-12": [ "ATLAS-SUSY-2013-05" ] }
        database.createLinksToCombinationsMatrix()
        self.assertTrue ( ds1.isCombinableWith ( ds3 ) )
        database.clearLinksToCombinationsMatrix()

    def testViaFields(self):
        """ check that combinability via "combinableWith" fields works 
        in the test base, ATLAS-SUSY-2013-05 is marked as combinable with ATLAS-SUSY-2013-02
        """
        er1 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-05" ] )[0]
        ds1 = er1.datasets[0]
        er2 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-02" ] )[0]
        ds2 = er2.datasets[0]
        er3 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-12" ] )[0]
        ds3 = er3.datasets[0]
        self.assertTrue ( ds1.isCombinableWith ( ds2 ) )
        self.assertTrue ( ds2.isCombinableWith ( ds1 ) )
        self.assertFalse ( ds1.isCombinableWith ( ds3 ) )
        self.assertTrue ( er1.isCombinableWith ( ds2 ) )

    def testAtDatasetLevel(self):
        """ check that things also work at the level of datasets
        """
        er1 = database.getExpResults( analysisIDs = [ "CMS-PAS-EXO-16-036" ] )[0]
        ds1 = er1.datasets[0]
        ds2 = er1.datasets[1]
        er2 = database.getExpResults( analysisIDs = [ "CMS-SUS-16-050-agg" ] )[0]
        ds3 = er2.datasets[0]
        ds4 = er2.datasets[1]
        self.assertFalse ( ds1.isCombinableWith ( ds3 ) )
        self.assertFalse ( ds1.isCombinableWith ( ds2 ) )
        self.assertFalse ( ds2.isCombinableWith ( ds3 ) )
        self.assertFalse ( ds1.isCombinableWith ( ds3 ) )
        self.assertFalse ( er1.isCombinableWith ( er2 ) )
        database.combinationsmatrix = { "CMS-PAS-EXO-16-036:c000": [ "CMS-SUS-16-050-agg:ar1" ] }
        database.createLinksToCombinationsMatrix()
        self.assertFalse ( ds1.isCombinableWith ( ds2 ) )
        self.assertFalse ( ds2.isCombinableWith ( ds3 ) )
        self.assertTrue ( ds1.isCombinableWith ( ds3 ) )
        self.assertFalse ( er1.isCombinableWith ( er2 ) )
        database.clearLinksToCombinationsMatrix()

        database.combinationsmatrix = { "CMS-SUS-16-050-agg:ar1": [ "CMS-PAS-EXO-16-036:c000", "CMS-PAS-EXO-16-036:c100" ] }
        database.createLinksToCombinationsMatrix()
        isc =  ds1.isCombinableWith ( ds2 )
        self.assertFalse ( ds1.isCombinableWith ( ds2 ) )
        self.assertTrue ( ds2.isCombinableWith ( ds3 ) )
        self.assertTrue ( ds1.isCombinableWith ( ds3 ) )
        self.assertFalse ( er1.isCombinableWith ( er2 ) )
        database.clearLinksToCombinationsMatrix()

        database.combinationsmatrix = { "CMS-PAS-EXO-16-036": [ "CMS-SUS-16-050-agg" ] }
        database.createLinksToCombinationsMatrix()
        self.assertFalse ( ds1.isCombinableWith ( ds2 ) )
        self.assertTrue ( ds2.isCombinableWith ( ds3 ) )
        self.assertTrue ( ds1.isCombinableWith ( ds3 ) )
        self.assertTrue ( ds1.isCombinableWith ( ds4 ) )
        self.assertTrue ( er1.isCombinableWith ( er2 ) )
        database.clearLinksToCombinationsMatrix()
        database.combinationsmatrix = { "CMS-PAS-EXO-16-036:c000": [ "CMS-PAS-EXO-16-036:c100" ] }
        database.createLinksToCombinationsMatrix()
        self.assertTrue ( ds1.isCombinableWith ( ds2 ) )
        database.clearLinksToCombinationsMatrix()

if __name__ == "__main__":
    unittest.main()
