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
from smodels.theory.particle import Particle
from smodels.theory.vertex import Vertex
from smodels.theory.branch import Branch,createBranchFromStr,decayBranches
from smodels.tools.physicsUnits import GeV, fb
import pickle

#Load the particle dictionaries
f = open("particleDefinitions.pcl","rb")
modelParticles = pickle.load(f)
particlesDict = dict([[p._name,p] for p in modelParticles])
f.close()


g = Particle(mass=500.*GeV,_pid=1000021, zParity=-1, _width = 1*GeV)
sq1 = Particle(mass = 150.*GeV, zParity = -1,_pid=1000005, _width = 1*GeV)
sq2 = Particle(mass = 100.*GeV, zParity = -1,_pid=2000005, _width = 1*GeV)
sq2b = Particle(mass = 100.*GeV, zParity = -1, _pid = 1000004, _width = 1*GeV)
sn1 = Particle(_name='neutralino1', mass = 50.*GeV, zParity = -1, _pid = 1000022, _width = 0.*GeV)
sn1B = Particle(_name='neutralino1', zParity = -1, _pid = 1000022, eCharge = 0, qColor = 0, _width = 0.*GeV)
u = Particle(_name='u', zParity = 1, _pid = 3, _width = 0.*GeV, eCharge = 2./3.)
d = Particle(_name='d',mass = 0.0*GeV, zParity = 1, _pid = 2, _width = 0.*GeV, eCharge = -1./3.)



class BranchTest(unittest.TestCase):
        
    def testBranch(self):
        
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v0B = Vertex(inParticle=None, outParticles=[sq1])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,d])
        v2b = Vertex(inParticle=sq1, outParticles=[sq2b,d])
        v3 = Vertex(inParticle=sq2, outParticles=[d,u,sn1])
        v4 = Vertex(inParticle=sq1, outParticles=[d,u,sn1])
        b1 = Branch(vertices = [v0,v1,v2,v3])
        b1b = Branch(vertices = [v0,v1,v2b,v3])
        b2 = Branch(vertices = [v0B,v2,v3])
        b3 = Branch(vertices = [v0,v1,v3])
        b4 = Branch(vertices = [v0,v1,v4])
        
        self.assertEqual( str(b1) == '[[d,u],[d],[d,u]]', True)
        self.assertEqual( b1 == b1b, True)
        self.assertEqual(b1 > b2, True)  #Larger by number of vertices
        self.assertEqual(b2 > b3, False)  #Smaller by number of outgoing particles
        self.assertEqual(b3 > b4, False)  #Smaller by mass of sq2
        self.assertEqual(len(b1) == 4, True)
        
    def testBranchInclusive(self):

        em = particlesDict['e-']
        mup = particlesDict['mu+']
        L = particlesDict['L']
        e = particlesDict['e']
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,em,em,mup])

        v1b = Vertex(inParticle=g, outParticles=[d,u,sq1])
        v2b = Vertex(inParticle=sq1, outParticles=[L,e,e,sq2])
        b1 = Branch(vertices = [v0,v1,v2])
        b1b = Branch(vertices = [v0,v1b,v2b])
        
        self.assertEqual( b1 == b1b, True) #Test if inclusive label comparison works
        
    def testBranchStr(self):
        
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,d])
        v3 = Vertex(inParticle=sq2, outParticles=[d,u,sn1B])
        b1 = Branch(vertices = [v0,v1,v2,v3])
                        
        bstr = createBranchFromStr('[[u,d],[d],[u,d]]',particlesDict)
        
        self.assertEqual( b1 == bstr, True)
        
    def testBranchAddVertex(self):
        
        v0 = Vertex(inParticle=None, outParticles=[g])
        v1 = Vertex(inParticle=g, outParticles=[sq1,u,d])
        v2 = Vertex(inParticle=sq1, outParticles=[sq2,d])
        v3 = Vertex(inParticle=sq2, outParticles=[d,u,sn1B])
        b1 = Branch(vertices = [v0,v1,v2])
        b2 = None
        try:
            b2 = b1._addVertex(v2)
        except SModelSTheoryError:
            pass
        self.assertEqual(b2 is None, True) #Should fail (incoming and outgoing do not match)
        
        b2 = b1._addVertex(v3)
        
        bstr = createBranchFromStr('[[u,d],[d],[u,d]]',particlesDict)
        self.assertEqual( b2 == bstr, True)
        

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
