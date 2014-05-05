#!/usr/bin/env python

"""
.. module:: testXSecComputer
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import logging
import setPath
from smodels.tools import xsecComputer
import logging.config

class XSecTest(unittest.TestCase):
    logging.config.fileConfig( "./logging.conf" )
    logger = logging.getLogger(__name__)

    def testLO(self):
        self.logger.info ( "test LO xsecs @ 8 TeV" )
        w = xsecComputer.computeXSec(8,0,1000, "../inputFiles/slha/andrePT4.slha").getDictionary()
        w8lo= 1000 * w[(1000023, 1000024)]['8 TeV (LO)'].asNumber() 
        self.assertAlmostEqual(w8lo, 35.014621117 )  ## 35.01 fb

    def testNLL (self):
        self.logger.info ( "test NLL xsecs @ 8 TeV" )
        filename="../inputFiles/slha/squarks.slha"
        w = xsecComputer.computeXSec(8,2,1000, filename ).getDictionary()
        w8nll = 1000 * w[(1000001, 1000002)]['8 TeV (NLO+NLL)'].asNumber() 
        self.assertAlmostEqual( w8nll, 60.915027554653705 ) ## 61 fb
if __name__ == "__main__":
    unittest.main()
