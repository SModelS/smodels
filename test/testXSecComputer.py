#!/usr/bin/env python

"""
.. module:: testXSecComputer
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

from smodels.tools import xsecComputer
from smodels.tools.xsecComputer import LO, NLL
from smodels.tools.physicsUnits import TeV, fb
import unittest
import logging
import logging.config

class XSecTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    logger = logging.getLogger(__name__)

    def testLOGlu(self):
        """ test the computation of LO cross section """
        self.logger.info ( "test LO xsecs @ 8 TeV" )
        slhafile="../inputFiles/slha/simplyGluino.slha"
        w = xsecComputer.computeXSec(8*TeV,LO,1000, slhafile ).getDictionary()
        w8lo= w[(1000021, 1000021)]['8 TeV (LO)'].asNumber( fb )
        self.assertAlmostEqual(w8lo, 268.4799000000022  ) ## 268.48 fb

    def testNLLGlu(self):
        """ test the computation of LO cross section """
        self.logger.info ( "test LO xsecs @ 8 TeV" )
        slhafile="../inputFiles/slha/simplyGluino.slha"
        w = xsecComputer.computeXSec(8*TeV,NLL,1000, slhafile ).getDictionary()
        w8lo= w[(1000021, 1000021)]['8 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual(w8lo, 583.1651907900048 ) ## 583.165 fb

    def testNLLSq (self):
        """ test the computation of NLL cross section """
        self.logger.info ( "test NLL xsecs @ 8 TeV" )
        filename="../inputFiles/slha/simplyGluino.slha"
        w = xsecComputer.computeXSec(8*TeV,NLL,1000, filename ).getDictionary()
        w8nll = w[(1000021, 1000021)]['8 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual( w8nll, 583.16519079 ) ## 10.95 pb


if __name__ == "__main__":
    unittest.main()
