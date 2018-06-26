#!/usr/bin/env python

"""
.. module:: testUpdateParticles
   :synopsis: Tests the function of smodels.theory.updateParticles 
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../../")
import unittest
from smodels.particleDefinitions import SMpdgs, BSMpdgs
from smodels.theory.updateParticles import getMassWidthBranches, addLongLived, calculateProbabilities, addPromptAndDisplaced
from smodels.theory.particleClass import Particles
from smodels.theory.branch import Branch
from smodels.tools.physicsUnits import eV, GeV
from smodels.theory.element import Element

n1 = Particles(Z2parity='odd', label='N1', pdg=1000022, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, branches=None) 
st1 = Particles(Z2parity='odd', label='st_1', pdg=1000006, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, branches=None)
gluino = Particles(Z2parity='odd', label='gluino', pdg=1000021, mass=None, eCharge=0, colordim=8, spin=1./2, width=None, branches=None)

slhafile = slhafile = '../../inputFiles/slha/simplyGluino.slha'
BSMList = [n1, st1, gluino]

getMassWidthBranches(slhafile, BSMList)

#METbranch = Branch()

class UpdateParticlesTest(unittest.TestCase):     
   
    def test1getMassWidthBranches(self):        
        
        self.assertEqual( n1.mass, 200.*GeV)    
        self.assertEqual( st1.width, 7.27235497E+00*GeV)
        self.assertEqual( gluino.branches[0].br, 0.5 )  
        self.assertEqual( gluino.branches[0].ids, [1000022, -1, 1] )  
        
        
    def test2addLongLived(self):
        
        gluino.width = 10**(-35)*GeV        
        addLongLived(BSMList)
        
        self.assertEqual( gluino.branches[0].br, 0. )
        self.assertEqual( gluino.branches[2].br, 1. ) 


    def test3calculateProbabilities(self):
      
        F_long, F_prompt, F_displaced = calculateProbabilities(gluino)
        self.assertEqual(F_long, 1.)
        self.assertEqual(F_prompt, 0.)
        self.assertEqual(F_displaced, 0.)
    
    
    def test4addPromptAndDisplaced(self):        
    
        branch1 = Branch()
        branch1.BSMparticles = [[n1]]
        probabilities1, branches1 = addPromptAndDisplaced(branch1)
        self.assertEqual(len(probabilities1), 1)
        self.assertEqual(probabilities1[0], 1.)
        self.assertEqual(branches1[0].decayType, 'METonly')
        
        branch2 = Branch()
        branch2.BSMparticles = [[gluino]]
        probabilities2, branches2 = addPromptAndDisplaced(branch2)
        self.assertEqual(len(probabilities2), 1)
        self.assertEqual(probabilities2[0], 1.)
        self.assertEqual(branches2[0].decayType, 'longlived')
        
        branch3 = Branch()
        branch3.BSMparticles = [[st1]]
        probabilities3, branches3 = addPromptAndDisplaced(branch3)
        self.assertEqual(len(probabilities3), 2)
        self.assertEqual(probabilities3[0], 1.)
        self.assertEqual(probabilities3[1], 0.)
        self.assertEqual(branches3[0].decayType, 'prompt') 
        self.assertEqual(branches3[1].decayType, 'displaced(neither jet nor lepton)')                
        

               
        
if __name__ == "__main__":
    unittest.main()        
