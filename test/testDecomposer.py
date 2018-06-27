#!/usr/bin/env python

"""
.. module:: testAsciiGraph
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.particleDefinitions import BSM
from smodels.theory.model import Model
from smodels.installation import installDirectory
from smodels.theory import decomposer
from smodels.theory.element import Element 
from smodels.tools.physicsUnits import GeV,pb,TeV

class DecomposerTest(unittest.TestCase):

    def testDecomposerLHE(self):

        filename = "%sinputFiles/lhe/simplyGluino.lhe" %(installDirectory())  
        model = Model(filename, BSM)
        model.updateParticles()
        
        topList = decomposer.decompose(model)
        self.assertTrue(len(topList.getElements()) == 1)
        element = topList.getElements()[0]
        el = Element("[[[q,q~]],[[q,q~]]]",finalState=['MET','MET'])
        self.assertTrue(el == element)
        bsmLabels = [[bsm.label for bsm in branch] for branch in element.getBSMparticles()]
        self.assertEqual(bsmLabels,[['gluino','N1']]*2)
        self.assertAlmostEqual(element.getMasses(),[[675.*GeV,200.*GeV]]*2)
        xsec = [xsec for xsec in element.weight if xsec.info.sqrts == 8.*TeV][0]
        xsec = xsec.value.asNumber(pb)        
        self.assertAlmostEqual(xsec,0.262,3)


    def testDecomposerSLHA(self):

        filename = "%sinputFiles/slha/simplyGluino.slha" %(installDirectory())  
        model = Model(filename, BSM)
        model.updateParticles()
        
        topList = decomposer.decompose(model)
        self.assertTrue(len(topList.getElements()) == 1)
        element = topList.getElements()[0]
        el = Element("[[[q,q~]],[[q,q~]]]",finalState=['MET','MET'])
        self.assertTrue(el == element)
        bsmLabels = [[bsm.label for bsm in branch] for branch in element.getBSMparticles()]
        self.assertEqual(bsmLabels,[['gluino','N1']]*2)
        self.assertAlmostEqual(element.getMasses(),[[675.*GeV,200.*GeV]]*2)
        xsec = [xsec for xsec in element.weight if xsec.info.sqrts == 8.*TeV][0]
        xsec = xsec.value.asNumber(pb)
        self.assertAlmostEqual(element.weight[0].value.asNumber(pb),0.572,3)



if __name__ == "__main__":
    unittest.main()