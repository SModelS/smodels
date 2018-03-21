#!/usr/bin/env python

"""
.. module:: testParticleClass
   :synopsis: Tests the smodels.theory.particleClass.Particles and .ParticleList class
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""
import sys
sys.path.insert(0,"../")
import unittest
import copy
from smodels.theory.particleClass import Particles, ParticleList, particleInList
from smodels.tools.physicsUnits import GeV

p1 = Particles(Z2parity='odd', label='p1', pdg=None, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p2 = Particles(Z2parity='odd', label='p1', pdg=1000021, mass=None, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p3 = Particles(Z2parity='odd', label='p3', pdg=1, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p4 = Particles(Z2parity='odd', label='p4', pdg=None, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)

p1c = copy.copy(p1)




class ParticleTest(unittest.TestCase):
    def testParticle(self):     
        #only comparisons by label for now           
        self.assertEqual( p1.label, 'p1')
        self.assertEqual( p1 == p2, True) 
        self.assertEqual( p1 > p3, False)     
        self.assertEqual( p4 == p3 == p2, False)
        self.assertEqual( p1c == p1, True)
        

    def testParticleList(self):
        l1 = ParticleList(label='plist', particles=[p1,p2])
        from smodels.SMparticleDefinitions import lList, LList

        self.assertEqual( l1.label == 'plist', True)
        self.assertEqual( lList > LList, False) #Llist is longer
        self.assertEqual( l1 == lList, False) 
        self.assertEqual( particleInList(p1,[l1]), True) #l1 contains p1
        
          
if __name__ == "__main__":
    unittest.main()
