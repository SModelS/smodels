#!/usr/bin/env python

"""
.. module:: testBranchClass
   :synopsis: Tests the theory.branch.Branch class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from smodels.theory.exceptions import SModelSTheoryError
import unittest
from smodels.theory.particle import Particle, ParticleList
from smodels.theory.branch import Branch, decayBranches, BranchWildcard
from smodels.tools.physicsUnits import GeV, fb, MeV
import pickle

#Load the particle dictionaries

g = Particle(mass=500.*GeV,pdg=1000021, Z2Parity=-1, width = 1*GeV)
sq1 = Particle(mass = 150.*GeV, Z2Parity = -1,pdg=1000005, width = 1*GeV)
sq2 = Particle(label = 'sq2', mass = 100.*GeV, Z2Parity = -1,pdg=2000005, colordim = 3, width = 1*GeV)
sq2b = Particle(label = 'sq2b', mass = 100.*GeV, Z2Parity = -1, pdg = 1000004, colordim = 3,  width = 1*GeV)
sn1 = Particle(label='neutralino1', mass = 50.*GeV, Z2Parity = -1, pdg = 1000022, width = 0.*GeV)
sn1B = Particle(label='neutralino1', Z2Parity = -1, pdg = 1000022, eCharge = 0, colordim = 0, width = 0.*GeV)
u = Particle(label='q', Z2Parity = 1, pdg = 3, width = 0.*GeV, eCharge = 2./3.)
e = Particle(Z2parity='even', label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, width = 0.*GeV, decays=[]) 



class BranchTest(unittest.TestCase):
        
    def testBranch(self):
        
        b1 = Branch(info = '[[q,q],[q],[q,q]]', finalState = 'MET')
        b2 = Branch(info = '[[q]]')
        b2.BSMparticles = [g]
        b3 = Branch(info = '[[e-]]')       
        b4 = Branch(info = '[[q]]')
        b4.BSMparticles = [sq1]
        b5 = BranchWildcard()

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
        b6.BSMparticles = [sq2]
        b6b = Branch( '[[e-],[q]]' )
        b6b.BSMparticles = [sq2b]
        
        b6.combineWith(b6b)
        
        self.assertTrue( isinstance(b6.BSMparticles[0],ParticleList) )
        self.assertEqual( b6.BSMparticles[0].label, 'BSM (combined)' )
        self.assertEqual(len(b6.BSMparticles[0].particles), 2)
        
        
        

    def testBranchDecay(self):
        
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v1b = Vertex(inParticle=g, outParticles=[sq2,u,d])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,d])
        v2b = Vertex(inParticle=sq1, outParticles=[sq2,u])            
        v3 = Vertex(inParticle=sq2, outParticles=[d,u,sn1])
        
        g._decayVertices = [v1,v1b]
        sq1._decayVertices = [v2,v2b]
        sq2._decayVertices = [v3]
        
        b0 = Branch(vertices = [v0])
        b1L = [Branch(vertices = [v0,v1]),Branch(vertices = [v0,v1b])]
        b2L = [Branch(vertices = [v0,v1,v2]),Branch(vertices = [v0,v1,v2b])]       
        
        d1L = b0.decay()
        self.assertEqual(len(d1L) == 2, True)
        self.assertEqual(sorted(d1L) == sorted(b1L), True)
        
        d2L = b1L[0].decay()
        self.assertEqual(len(d2L) == 2, True)
        self.assertEqual( d2L == b2L, True)                     
        
        
        
    def testDecayBranches(self):
        
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v1b = Vertex(inParticle=g, outParticles=[sq2,u,d])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,d])
        v2b = Vertex(inParticle=sq1, outParticles=[sq2,u])            
        v3 = Vertex(inParticle=sq2, outParticles=[d,u,sn1])
        
        v1.br = 0.8
        v1b.br = 0.5
        v2.br = 0.1
        v2b.br = 0.9
        v3.br = 0.3
        
        g._decayVertices = [v1,v1b]
        sq1._decayVertices = [v2,v2b]
        sq2._decayVertices = [v3]
        
        
        b0 =  Branch(vertices = [v0])
        b0.maxWeight = 10.*fb
        dList = decayBranches([b0])
        self.assertEqual(len(dList) == 3, True)
        
        bList = [Branch(vertices = [v0,v1,v2,v3]),Branch(vertices = [v0,v1,v2b,v3])
                 ,Branch(vertices = [v0,v1b,v3])]
        bList[0].maxWeight = b0.maxWeight*v1.br*v2.br*v3.br
        bList[1].maxWeight = b0.maxWeight*v1.br*v2b.br*v3.br
        bList[2].maxWeight = b0.maxWeight*v1b.br*v3.br             
        
        
        self.assertEqual(sorted(bList) == sorted(dList), True)
        
        weightsA = [b.maxWeight for b in sorted(bList)]
        weightsB = [b.maxWeight for b in sorted(dList)]
        self.assertEqual(weightsA == weightsB, True)
        
        
if __name__ == "__main__":
    unittest.main()
