#!/usr/bin/env python

"""
.. module:: testUpperLimit
   :synopsis: Test smsInterpolation.upperLimit with various inputs

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
from smodels.tools.physicsUnits import GeV, pb
from smodels.experiment.smsInterpolation import upperLimit
from smodels.experiment import smsHelpers
from smodels.installation import installDirectory

## smsHelpers.base = installDirectory() + 'test/database/'
smsHelpers.base =  './database/'

class UpperLimitTest(unittest.TestCase):

    def testDirectDecay(self):
        ul = float(upperLimit("SUS12028","T2",[400*GeV,100*GeV])/pb)
        self.assertAlmostEquals ( ul, 0.43757864832878113)

    def testCascadeDecay(self):
        ul = float(upperLimit("ATLAS_CONF_2013_048","T6bbWW",[500*GeV,400*GeV,100*GeV])/pb)
        self.assertAlmostEquals( ul , 0.10152592613065325 )

    def testMissingUnits(self):
        ul = upperLimit("ATLAS_CONF_2013_048", "T6bbWW", [500, 400, 100])
        self.assertTrue(ul==None)

if __name__ == "__main__":
    unittest.main()

