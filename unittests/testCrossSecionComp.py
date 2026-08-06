#!/usr/bin/env python3

"""
.. module:: testCrossSectionComp
   :synopsis: Tests xsec comparison.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.base.crossSection import XSection,XSectionInfo,XSectionList
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, pb


class XSecCompTest(unittest.TestCase):

    def testCrossSectionComp(self):
        
        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,
                                    ignorePromptQNumbers=['spin','eCharge','colordim'])

        xsecA = model.xsections[0]
        xsecB = model.xsections[1]
        self.assertEqual(xsecA.info.label,'8 TeV (LO)')
        self.assertEqual(xsecB.info.label,'13 TeV (LO)')

        self.assertEqual(xsecA.pid,(-1000037, 1000023))
        self.assertEqual(xsecB.pid,(-1000037, 1000023))

        self.assertTrue(xsecA < xsecB)
        self.assertFalse(xsecA > 5*fb)
        self.assertTrue(xsecA > 4e-3)

        maxXsec = model.xsections.getMaxXsec()
        self.assertAlmostEqual(maxXsec.asNumber(pb),11.3,1)
        
        self.assertFalse(model.xsections > 15)
        self.assertTrue(model.xsections > 15.*fb)
        self.assertFalse(model.xsections < xsecB)
        
if __name__ == "__main__":
    unittest.main()
