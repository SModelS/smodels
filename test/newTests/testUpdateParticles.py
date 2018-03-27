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
from smodels.theory.updateParticles import getMassWidthBranches, getPromptDecays
from smodels.theory.particleClass import Particles
from smodels.tools.physicsUnits import eV, GeV

n1 = Particles(Z2parity='odd', label='N1', pdg=1000022, mass=None, eCharge=0, colordim=0, spin=1./2, width=None, branches=None) 
st1 = Particles(Z2parity='odd', label='st_1', pdg=1000006, mass=None, eCharge=-1./3, colordim=3, spin=0, width=None, branches=None)
gluino = Particles(Z2parity='odd', label='gluino', pdg=1000021, mass=None, eCharge=0, colordim=8, spin=1./2, width=None, branches=None)

slhafile = slhafile = '../../inputFiles/slha/simplyGluino.slha'
BSMList = [n1, st1, gluino]

getMassWidthBranches(slhafile, BSMList)

class UpdateParticlesTest(unittest.TestCase):     
   
    def testgetMassWidthBranches(self):        
        
        self.assertEqual( n1.mass, 200.*GeV)    
        self.assertEqual( st1.width, round(7.27235497, 1)*GeV)
        self.assertEqual( gluino.branches[0].br, 0.5 )  
        self.assertEqual( gluino.branches[0].ids, [1000022, -1, 1] )  
        
        
    def testgetPromptDecays(self):

        gluino.width = 10**(-17)*eV        
        getPromptDecays(BSMList)
        
        self.assertEqual( round(gluino.branches[0].br, 5), 0. )
        self.assertEqual( round(gluino.branches[2].br, 5), 1. ) 
        
        # if the width is smaller than ~10^-2 eV, the BR will be reweighted with F_prompt (but no stable branch will be added)        
        st1.width = 10**(-2.2)*eV
        getPromptDecays(BSMList)
        self.assertEqual( round(st1.branches[0].br, 5), 1. )

        # if the width is smaller than ~10^-5 eV, also a stable branch will be added          
        st1.width = 10**(-5.0545)*eV
        getPromptDecays(BSMList)
        self.assertEqual( round(st1.branches[1].br), 0. )          
               
        
if __name__ == "__main__":
    unittest.main()        
