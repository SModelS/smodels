#!/usr/bin/env python3

"""
.. module:: testReweighting
   :synopsis: Tests the function of lifetime reweighting
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models import SMparticles, mssm
from smodels.theory.branch import Branch
from smodels.tools.reweighting import calculateProbabilities, addPromptAndDisplaced
from smodels.tools.physicsUnits import GeV

class ReweightingTest(unittest.TestCase):     
   
    def testcalculateProbabilities(self):

        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1.*10**(-30)*GeV
        prob = calculateProbabilities(gluino.totalwidth.asNumber(GeV))
        F_long, F_prompt, F_displaced = prob['F_long'],prob['F_prompt'],prob['F_displaced']
        self.assertAlmostEqual(F_long, 1.)
        self.assertEqual(F_prompt, 0.)
        self.assertAlmostEqual(F_displaced, 0.)
    
    
    def testaddPromptAndDisplaced(self):

        n1 = mssm.n1.copy()
        n1.totalwidth = 0.*GeV
        st1 = mssm.st1.copy()
        st1.totalwidth = 2.*GeV
        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1.*10**(-30)*GeV
        t = SMparticles.t
    
        branch1 = Branch()
        branch1.oddParticles = [n1]
        probabilities1, branches1 = addPromptAndDisplaced(branch1)
        self.assertEqual(len(probabilities1), 1)
        self.assertEqual(probabilities1[0], 1.)
        self.assertEqual(branches1[0]._decayType, 'METonly')
        
        branch2 = Branch()
        branch2.oddParticles = [gluino]
        probabilities2, branches2 = addPromptAndDisplaced(branch2)
        self.assertEqual(len(probabilities2), 1)
        self.assertAlmostEqual(probabilities2[0], 1.)
        self.assertEqual(branches2[0]._decayType, 'longlived')
        
        branch3 = Branch()
        branch3.oddParticles = [st1,n1]
        branch3.evenParticles = [[t]]
        probabilities3, branches3 = addPromptAndDisplaced(branch3)
        self.assertEqual(len(probabilities3), 2)
        self.assertEqual(probabilities3[0], 1.)
        self.assertEqual(probabilities3[1], 0.)
        self.assertEqual(branches3[0]._decayType, 'prompt') 
        self.assertEqual(branches3[1]._decayType, 'displaced(neither jet nor lepton)')                
        

               
        
if __name__ == "__main__":
    unittest.main()        
