#!/usr/bin/env python

"""
.. module:: testParticleComparison
   :synopsis: Tests the functions of smodels.theory.particleComparison 
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../../")
import unittest
from smodels.theory.particleClass import Particles
from smodels.theory.particleComparison import compareBSMparticles, simParticles, ParticleWildcard
from smodels.tools.physicsUnits import GeV

p1 = Particles(Z2parity='odd', label='p1', pdg=None, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p2 = Particles(Z2parity='odd', label='p2', pdg=1000021, mass=100.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p3 = Particles(Z2parity='odd', label='p3', pdg=1, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p4 = Particles(Z2parity='odd', label='p4', pdg=None, mass=110.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)

p5 = Particles(Z2parity='odd', label='p5', pdg=1, mass=50.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)
p6 = Particles(Z2parity='odd', label='p6', pdg=None, mass=60.*GeV, eCharge=None, colordim=None, spin=None, width=None, branches=None)



class ParticleComparisonTest(unittest.TestCase):
   
    def testBSMcomparison(self):     
        BSM1 = [[p1, p3],[p2, p4]]
        BSM2 = [[p5,p6]]
        m1m2eq, comp = compareBSMparticles(BSM1, BSM2)

        self.assertEqual( m1m2eq, False)
        self.assertEqual( comp, True)
        
        
    def testSimParticles(self):
        from smodels.SMparticleDefinitions import lList, e, mu
        l1 = [e,mu]
        l2 = [lList, lList]
        
        self.assertEqual( simParticles(l1, l2), True)
        
    def testParticleWildCard(self):
        particle2 = ParticleWildcard()
        
        self.assertTrue(isinstance(p1, Particles))
        self.assertTrue(p1 == particle2)
        self.assertTrue(particle2 == p1)
        
        
        
        
        
if __name__ == "__main__":
    unittest.main()        
