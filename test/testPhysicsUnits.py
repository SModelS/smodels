#!/usr/bin/env python

"""
.. module:: testPhysicsUnits
   :synopsis: Tests conversion of the physics units.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
from smodels.tools.physicsUnits import addunit, rmvunit, GeV, TeV, eV, fb

class PhysicsUnitsTest(unittest.TestCase):
    def testXSecs(self):
        x1=addunit ( 3000.0, "fb" )
        x2=addunit ( 3.0, "pb" )
        self.assertAlmostEqual(x1,x2)
    def testEnergies(self):
        d1= 7000.0 * GeV
        d2= 7.0 * TeV
        self.assertAlmostEqual(d1.asNumber(GeV),d2.asNumber(GeV))

    def testLumi(self):
        d1=addunit ( 2.0, "fb-1" )
        self.assertAlmostEqual ( (d1*.5*fb).asNumber(), 1.0 )


if __name__ == "__main__":
    unittest.main()
