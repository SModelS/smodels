#!/usr/bin/env python

"""
.. module:: NllFast test
   :synopsis: Tests nll computations via nllfast.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
import math
from smodels.tools import nllFast, externalNllFast
from smodels.tools.physicsUnits import TeV

class NllFastTest(unittest.TestCase):
    def testGluino7 (self):
        o = nllFast.getKfactorsFor ( (1000021, 1000021 ), 7*TeV, "../inputFiles/slha/simplyGluino.slha" )
        self.assertAlmostEqual ( o[0], 2.1 )
        self.assertAlmostEqual ( o[1], 1.08 )
    def testSquark8 (self):
        o = nllFast.getKfactorsFor ( (1000001, 1000001 ), 8*TeV, "../inputFiles/slha/gluino_squarks.slha" )
        self.assertAlmostEqual ( o[0], 1.22 )
        self.assertAlmostEqual ( o[1], 1.02 )
    def testSquark13 (self):
        o = nllFast.getKfactorsFor ( (1000001, 1000001 ), 13*TeV, "../inputFiles/slha/gluino_squarks.slha" )
        self.assertAlmostEqual ( o[0], 1.24 )
        self.assertAlmostEqual ( o[1], 1.01 )
    def testWeakino8 (self):
        o = nllFast.getKfactorsFor ( (1000022, 1000022 ), 13*TeV, "../inputFiles/slha/complicated.slha" )
        self.assertEqual ( o[0], None )
        self.assertEqual ( o[1], None )

if __name__ == "__main__":
    unittest.main()
