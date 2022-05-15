#!/usr/bin/env python3

"""
.. module:: testPhysicsUnits
   :synopsis: Tests conversion of the physics units.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.physicsUnits import GeV, TeV, fb, pb

class PhysicsUnitsTest(unittest.TestCase):
    def testXSecs(self):
        x1= 3000.0 * fb
        x2= 3.0 * pb
        self.assertAlmostEqual( x1.asNumber(fb),x2.asNumber(fb) )
    def testEnergies(self):
        d1= 7000.0 * GeV
        d2= 7.0 * TeV
        self.assertAlmostEqual(d1.asNumber(GeV),d2.asNumber(GeV))

    def testLumi(self):
        d1= 2.0 / fb 
        self.assertAlmostEqual ( (d1*.5*fb).asNumber(), 1.0 )


if __name__ == "__main__":
    unittest.main()
