#!/usr/bin/env python

"""
.. module:: testReweighting
   :synopsis: Tests the function of smodels.theory.updateParticles 
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.particle import Particle
from smodels.share.models import SMparticles, mssm
from smodels.theory.branch import Branch
from smodels.theory.reweighting import calculateProbabilities, addPromptAndDisplaced
from smodels.tools.physicsUnits import eV, GeV
from smodels.theory.element import Element

n1 = mssm.n1
n1.totalwidth = 0.*GeV
st1 = mssm.st1
st1.totalwidth = 2.*GeV 
gluino = mssm.gluino
gluino.totalwidth = 1.*10**(-30)*GeV
t = SMparticles.t

class ReweightingTest(unittest.TestCase):     
   
    def testcalculateProbabilities(self):
      
        F_long, F_prompt, F_displaced = calculateProbabilities(gluino)
        self.assertAlmostEqual(F_long, 1.)
        self.assertEqual(F_prompt, 0.)
        self.assertAlmostEqual(F_displaced, 0.)
    
    
    def testaddPromptAndDisplaced(self):        
    
        branch1 = Branch()
        branch1.BSMparticles = [n1]
        probabilities1, branches1 = addPromptAndDisplaced(branch1)
        self.assertEqual(len(probabilities1), 1)
        self.assertEqual(probabilities1[0], 1.)
        self.assertEqual(branches1[0]._decayType, 'METonly')
        
        branch2 = Branch()
        branch2.BSMparticles = [gluino]
        probabilities2, branches2 = addPromptAndDisplaced(branch2)
        self.assertEqual(len(probabilities2), 1)
        self.assertAlmostEqual(probabilities2[0], 1.)
        self.assertEqual(branches2[0]._decayType, 'longlived')
        
        branch3 = Branch()
        branch3.BSMparticles = [st1,n1]
        branch3.particles = [[t]]
        probabilities3, branches3 = addPromptAndDisplaced(branch3)
        self.assertEqual(len(probabilities3), 2)
        self.assertEqual(probabilities3[0], 1.)
        self.assertEqual(probabilities3[1], 0.)
        self.assertEqual(branches3[0]._decayType, 'prompt') 
        self.assertEqual(branches3[1]._decayType, 'displaced(neither jet nor lepton)')                
        

               
        
if __name__ == "__main__":
    unittest.main()        
