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
from smodels.share.models.MSSMparticles import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.installation import installDirectory
from smodels.theory import decomposer
from smodels.theory.element import Element 
from smodels.tools.physicsUnits import GeV,pb,TeV,fb

class DecomposerTest(unittest.TestCase):

    def testDecomposerLHE(self):
 
        filename = "%sinputFiles/lhe/simplyGluino.lhe" %(installDirectory())  
        model = Model(BSMList,SMList,filename)
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
        model = Model(BSMList,SMList,filename)
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

 
    def testDecomposerLongLived(self):
  
        filename = "%sinputFiles/slha/longLived.slha" %(installDirectory())
        #Consider a simpler model
        newModel = [ptc for ptc in BSMList if not isinstance(ptc.pdg,list) and abs(ptc.pdg) in [1000015,1000022]] 
        model = Model(newModel,SMList,filename)
        model.updateParticles()
          
        topList = decomposer.decompose(model)
        self.assertTrue(len(topList.getElements()) == 10)
        expectedWeights = {str(sorted([['N1'],['N1']])).replace(' ','') : 0.020,
                           str(sorted([['sta_1'],['sta_1~']])).replace(' ','') : 0.26,
                           str(sorted([['sta_1'],['sta_1~','N1~']])).replace(' ','') : 0.13,
                           str(sorted([['sta_1~'],['sta_1','N1']])).replace(' ','') : 0.13,
                           str(sorted([['sta_1~','N1~'],['sta_1','N1']])).replace(' ','') : 0.065}
            
        for el in topList.getElements():
            bsmLabels = str(sorted([[bsm.label for bsm in branch] for branch in el.getBSMparticles()]))
            bsmLabels = bsmLabels.replace(' ','')
            xsec = el.weight.getXsecsFor(8.*TeV)[0].value.asNumber(fb)
            self.assertAlmostEqual(expectedWeights[bsmLabels], xsec,2)



if __name__ == "__main__":
    unittest.main()