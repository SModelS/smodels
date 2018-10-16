#!/usr/bin/env python3

"""
.. module:: testBranchClass
   :synopsis: Tests the theory.branch.Branch class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.particle import Particle, MultiParticle
from smodels.theory.branch import Branch, decayBranches, InclusiveBranch
from smodels.tools.physicsUnits import GeV, fb, MeV
import pyslha

#Load the particle dictionaries
    
sq1 = Particle(Z2parity='odd', mass = 150.*GeV, pdg=1000005, totalwidth = 1*GeV)
sq2 = Particle(Z2parity='odd',label = 'sq2', mass = 100.*GeV, pdg=2000005, colordim = 3, totalwidth = 1*GeV)
sq2b = Particle(Z2parity='odd',label = 'sq2b', mass = 100.*GeV,  pdg = 1000004, colordim = 3,  totalwidth = 1*GeV)
sn1 = Particle(Z2parity='odd',label='neutralino1', mass = 50.*GeV,  pdg = 1000022, width = 0.*GeV, decays=[])
sn1B = Particle(Z2parity='odd', label='neutralino1', pdg = 1000022, eCharge = 0, colordim = 0, totalwidth = 0.*GeV)
u = Particle(Z2parity='even', label='q',  pdg = 3, totalwidth = 0.*GeV, eCharge = 2./3.)
e = Particle(Z2parity='even', label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[]) 

g_decays = [pyslha.Decay(br=0.3,nda=3,ids=[-1,1,1000022],parentid=1000021),
                      pyslha.Decay(br=0.7,nda=3,ids=[-2,2,1000022],parentid=1000021)]
g_decays[0].daughters = [u,u,sn1]
g_decays[1].daughters = [u,u,sn1]
g = Particle(mass=500.*GeV,pdg=1000021, Z2Parity=-1, totalwidth = 1*GeV, decays = g_decays )

class BranchTest(unittest.TestCase):
        
    def testBranch(self):
        
        b1 = Branch(info = '[[q,q],[q],[q,q]]', finalState = 'MET')
        b2 = Branch(info = '[[q]]')
        b2.oddParticles = [g]
        b3 = Branch(info = '[[e-]]')       
        b4 = Branch(info = '[[q]]')
        b4.oddParticles = [sq1]
        b5 = InclusiveBranch()

        self.assertTrue( str(b1) == '[[q,q],[q],[q,q]]')
        self.assertTrue(b1 > b2)  #Larger by number of vertices
        self.assertEqual(b1, b5) #always true because b3 is a wildcard
        self.assertFalse(b2 < b3)  #Bigger by label of particles
        self.assertFalse(b2 < b4)  #Bigger by mass of g compared to sq1

        
    def testBranchInclusive(self):

        bi1 = Branch( '[[e-],[q]]' )
        bi1b = Branch( '[[l],[jet]]' )
        
        self.assertTrue( bi1.particlesMatch(bi1b)) #Test if inclusive label comparison works

        
    def testcombineWith(self):
        
        b6 = Branch( '[[e-],[q]]' )
        b6.oddParticles = [sq2]
        b6b = Branch( '[[e-],[q]]' )
        b6b.oddParticles = [sq2b]
        
        b6.combineWith(b6b)
        
        self.assertTrue( isinstance(b6.oddParticles[0],MultiParticle) )
        self.assertEqual( b6.oddParticles[0].label, 'BSM (combined)' )
        self.assertEqual(len(b6.oddParticles[0].particles), 2)
        
        
        

    def testBranchDecay(self):

        b = Branch(info = '[[q]]', finalState = 'MET')
        b.oddParticles = [g]
        
        newBranches = b.decayDaughter()
        self.assertEqual(len(newBranches), 2)
        self.assertEqual(newBranches[0].evenParticles, [[u],[u,u]])
        self.assertEqual(newBranches[0].oddParticles, [g,sn1])

        
        
        
        
    def testDecayBranches(self):
        
        b = Branch(info = '[[q]]', finalState = 'MET')
        b.maxWeight = 2.0*fb
        b.oddParticles = [g]
        
        dList = decayBranches([b])
        
        self.assertEqual(len(dList), 2)
        self.assertEqual(dList[0].maxWeight, .6*fb)
        self.assertEqual(dList[1].maxWeight, 1.4*fb)
        
        
if __name__ == "__main__":
    unittest.main()
