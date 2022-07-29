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
        ds1 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-05" ] )[0].datasets[0]
        ds2 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-02" ] )[0].datasets[0]
        ds3 = database.getExpResults( analysisIDs = [ "ATLAS-SUSY-2013-12" ] )[0].datasets[0]
        self.assertTrue ( ds1.isCombinableWith ( ds2 ) )
        self.assertTrue ( ds2.isCombinableWith ( ds1 ) )
        self.assertFalse ( ds1.isCombinableWith ( ds3 ) )


if __name__ == "__main__":
    unittest.main()
