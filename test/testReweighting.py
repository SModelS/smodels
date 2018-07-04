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
from smodels.theory.branch import Branch
from smodels.theory.reweighting import calculateProbabilities, addPromptAndDisplaced
from smodels.tools.physicsUnits import eV, GeV
from smodels.theory.element import Element

n1 = Particle(Z2parity='odd', label='N1', pdg=1000022, mass=None, eCharge=0, colordim=0, spin=1./2, totalwidth=0.*GeV, branches=None) 
st1 = Particle(Z2parity='odd', label='st_1', pdg=1000006, mass=None, eCharge=-1./3, colordim=3, spin=0, totalwidth=2.*GeV, branches=None)
gluino = Particle(Z2parity='odd', label='gluino', pdg=1000021, mass=None, eCharge=0, colordim=8, spin=1./2, totalwidth=1.*10**(-30)*GeV, branches=None)


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
        branch3.BSMparticles = [st1]
        probabilities3, branches3 = addPromptAndDisplaced(branch3)
        self.assertEqual(len(probabilities3), 2)
        self.assertEqual(probabilities3[0], 1.)
        self.assertEqual(probabilities3[1], 0.)
        self.assertEqual(branches3[0]._decayType, 'prompt') 
        self.assertEqual(branches3[1]._decayType, 'displaced(neither jet nor lepton)')                
        

               
        
if __name__ == "__main__":
    unittest.main()        
