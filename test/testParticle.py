#!/usr/bin/env python

"""
.. module:: testParticleClass
   :synopsis: Tests the smodels.theory.particleClass.Particles and .ParticleList class
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../")
import unittest
import copy
from smodels.theory.particle import Particle, ParticleList
from smodels.tools.physicsUnits import GeV
from smodels.theory.particleNames import getObjectFromLabel

p1 = Particle(Z2parity='odd', label='p1', pdg=None, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p2 = Particle(Z2parity='odd', label='p1', pdg=1000021, mass=None, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p3 = Particle(Z2parity='odd', label='p3', pdg=1, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p4 = Particle(Z2parity='odd', label='p4', pdg=None, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p5 = Particle(Z2parity='odd', label='p5+', pdg=2, mass=110.*GeV, eCharge=1., colordim=None, spin=None, width=None, branches=None)

p1c = copy.copy(p1)
p5cc = p5.chargeConjugate()
p3cc = p3.chargeConjugate()


class ParticleTest(unittest.TestCase):
    def testParticle(self):     
        #only comparisons by label for now           
        self.assertEqual( p1.label, 'p1')
        self.assertEqual( p1 , p2) 
        self.assertLess( p1 , p3)     
        self.assertFalse( p4 == p3 == p2) 
        self.assertEqual( p1c , p1)
        

    def testParticleList(self):
        l1 = ParticleList(label='plist', particles=[p1,p2])
        from smodels.SMparticleDefinitions import lList, LList

        self.assertEqual( l1.label, 'plist')    
        self.assertGreater(LList, lList)  #Llist is longer
        self.assertNotEqual( l1 , lList) 
        self.assertTrue(l1 == p1)
        self.assertTrue(sorted(l1.pdg) == sorted([1000021,None]))
        
    def testChargeConjugation(self):
        self.assertEqual(p5cc.label , 'p5-')
        self.assertEqual(p5cc.eCharge , -1.)
        self.assertEqual(p3cc.label, p3.label)
        self.assertEqual(p3cc.pdg, -1)
        
    def testParticleWildCard(self):
        particle2 = getObjectFromLabel('*')
        l1 = ParticleList(label='plist', particles=[p1,p2,p4])
        from smodels.SMparticleDefinitions import lList
        
        
        self.assertTrue(isinstance(p1, Particle))
        self.assertTrue(p1 == particle2)
        self.assertTrue(particle2 == p1)
        self.assertTrue(particle2 == l1)
        self.assertTrue(l1 == particle2)
        self.assertTrue(particle2 == lList)
        self.assertTrue(lList == particle2)

        
        
          
if __name__ == "__main__":
    unittest.main()
