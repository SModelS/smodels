#!/usr/bin/env python3

"""
.. module:: testBranchClass
   :synopsis: Tests the theory.branch.Branch class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.particle import Particle, MultiParticle,ParticleList
from smodels.theory.branch import Branch, decayBranches, InclusiveBranch
from smodels.tools.physicsUnits import GeV, MeV
from smodels.experiment.defaultFinalStates import finalStates
import pyslha

#Load the particle dictionaries

sq1 = Particle(Z2parity=-1, mass = 150.*GeV, pdg=1000005, totalwidth = 1*GeV)
sq2 = Particle(Z2parity=-1,label = 'sq2', mass = 100.*GeV, pdg=2000005, colordim = 3, totalwidth = 1*GeV)
sq2b = Particle(Z2parity=-1,label = 'sq2b', mass = 100.*GeV,  pdg = 1000004, colordim = 3,  totalwidth = 1*GeV)
sn1 = Particle(Z2parity=-1,label='neutralino1', mass = 50.*GeV,  pdg = 1000022, width = 0.*GeV, decays=[])
sn1B = Particle(Z2parity=-1, label='neutralino1', pdg = 1000022, eCharge = 0, colordim = 0, totalwidth = 0.*GeV)
u = Particle(Z2parity=1, label='q',  pdg = 3, totalwidth = 0.*GeV, eCharge = 2./3.)
e = Particle(Z2parity=1, label='e-', pdg=11, mass=0.5*MeV, eCharge=-1, colordim=0, spin=1./2, totalwidth = 0.*GeV, decays=[])

g_decays = [pyslha.Decay(br=0.3,nda=3,ids=[-1,1,1000022],parentid=1000021),
                      pyslha.Decay(br=0.7,nda=3,ids=[-2,2,1000022],parentid=1000021)]
g_decays[0].oddParticles = [sn1]
g_decays[0].evenParticles = ParticleList([u,u])
g_decays[1].oddParticles = [sn1]
g_decays[1].evenParticles = ParticleList([u,u])
g = Particle(mass=500.*GeV,pdg=1000021, Z2Parity=-1, totalwidth = 1*GeV, decays = g_decays )

class BranchTest(unittest.TestCase):

    def testBranch(self):

        b1 = Branch(info = '[[q,q],[q],[q,q]]', finalState = 'MET', model=finalStates)
        b2 = Branch(info = '[[q]]', model=finalStates)
        b2.oddParticles = [g]
        b3 = Branch(info = '[[e-]]', model=finalStates)
        b4 = Branch(info = '[[q]]', model=finalStates)
        b4.oddParticles = [sq1]
        b5 = InclusiveBranch(model=finalStates)

        self.assertTrue( str(b1) == '[[q,q],[q],[q,q]]')
        self.assertTrue(b1 > b2)  #Larger by number of vertices
        self.assertEqual(b1, b5) #always true because b3 is a wildcard
        self.assertTrue(b2 < b3)  #Bigger by odd particles
        self.assertFalse(b2 < b4)  #Bigger by mass of g compared to sq1


    def testBranchComp(self):

        gluino = Particle(mass=500.*GeV,pdg=1000021, Z2Parity=-1, eCharge=0, colordim=8, totalwidth = 1*GeV, decays = g_decays )
        sdR = Particle(mass=700.*GeV,pdg=2000001, Z2Parity=-1, eCharge=-1/3, colordim=3, totalwidth = 1*GeV, decays = None)
        sdL = Particle(mass=705.*GeV,pdg=1000001, Z2Parity=-1, eCharge=-1/3, colordim=3, totalwidth = 1*GeV, decays = None)
        ssL = Particle(mass=705.*GeV,pdg=1000003, Z2Parity=-1, eCharge=-1/3, colordim=3, totalwidth = 1*GeV, decays = None)
        suL = Particle(mass=705.*GeV,pdg=1000002, Z2Parity=-1, eCharge=2/3, colordim=3, totalwidth = 1*GeV, decays = None)
        C1m = Particle(mass=205.*GeV,pdg=1000003, Z2Parity=-1, eCharge=-1, colordim=1, totalwidth = 1*GeV, decays = None)
        N1 = Particle(mass=705.*GeV,pdg=1000022, Z2Parity=-1, eCharge=0, colordim=1)

        b1 = Branch('[[q], [b, t+], [W-]]', model=finalStates)
        b2 = Branch('[[q], [b, t+], [W-]]', model=finalStates)
        b3 = Branch('[[q], [b, t+], [W-]]', model=finalStates)
        b4 = Branch('[[q], [b, t+], [W-]]', model=finalStates)
        b1.oddParticles = [sdR,gluino,C1m,N1]
        b2.oddParticles = [sdL,gluino,C1m,N1]
        b3.oddParticles = [ssL,gluino,C1m,N1]
        b4.oddParticles = [suL,gluino,C1m,N1]

        self.assertNotEqual(b1,b2) #Same quantum numbers distinct masses
        self.assertEqual(b2,b3)  #Same masses and quantum numbers
        self.assertNotEqual(b3,b4) #Same masses distinct quantum number


    def testBranchInclusive(self):

        bi1 = Branch('[[e-],[q]]', model=finalStates)
        bi1b = Branch('[[l],[jet]]', model=finalStates)

        self.assertTrue(bi1 == bi1b) #Test if inclusive label comparison works


    def testcombineWith(self):

        b6 = Branch( '[[e-],[q]]', model=finalStates)
        b6.oddParticles = [sq2]
        b6b = Branch( '[[e-],[q]]', model=finalStates)
        b6b.oddParticles = [sq2b]

        b6 += b6b

        self.assertTrue(isinstance(b6.oddParticles[0],MultiParticle) )
        self.assertEqual(b6.oddParticles[0].label, 'sq2/sq2b' )
        self.assertEqual(len(b6.oddParticles[0].particles), 2)


    def testBranchDecay(self):

        b = Branch(info = '[[q]]', finalState = 'MET', model=finalStates)
        b.oddParticles = [g]

        newBranches = b.decayDaughter()
        self.assertEqual(len(newBranches), 2)
        self.assertEqual(newBranches[0].evenParticles, [ParticleList([u]),ParticleList([u,u])])
        self.assertEqual(newBranches[0].oddParticles, [g,sn1])


    def testDecayBranches(self):

        b = Branch(info = '[[q]]', finalState = 'MET', model=finalStates)
        b.maxWeight = 2.0
        b.oddParticles = [g]

        dList = decayBranches([b])

        self.assertEqual(len(dList), 2)
        self.assertEqual(dList[0].maxWeight, .6)
        self.assertEqual(dList[1].maxWeight, 1.4)


if __name__ == "__main__":
    unittest.main()
