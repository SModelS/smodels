#!/usr/bin/env python

"""
.. module:: testVertexClass
   :synopsis: Tests the theory.vertex.Vertex class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.particle import Particle
from smodels.theory.vertex import Vertex, createVertexFromStr
from smodels.tools.physicsUnits import GeV
import pickle

#Load the particle dictionaries
f = open("particleDefinitions.pcl","rb")
modelParticles = pickle.load(f)
particlesDict = dict([[p._name,p] for p in modelParticles])
f.close()


class VertexTest(unittest.TestCase):        
    def testVertex(self):
        
        p1 = Particle(mass=100.*GeV, zParity=-1)
        p3 = Particle( _pid=1, mass=110.*GeV, zParity=-1)
        p4 = Particle(mass = 110.*GeV, zParity = -1)
        
        sq1 = Particle(_name='squark1', mass = 110.*GeV, zParity = -1)
        sq2 = Particle(_name='squark1', mass = 100.*GeV, zParity = -1)
        u = Particle(_name='u', zParity = 1, eCharge = 2./3.)
        d = Particle(_name='d',mass = 0.01*GeV, zParity = 1, eCharge = -1./3.)
        
        v1 = Vertex(inParticle=p1, outParticles=[sq1,u,d])
        v2 = Vertex(inParticle=p1, outParticles=[sq1,d])
        v3 = Vertex(inParticle=p1, outParticles=[sq2,d])
        v5 = Vertex(inParticle=p3, outParticles=[sq1,u,d])
        v5b = Vertex(inParticle=p4, outParticles=[sq1,u,d])
        
        self.assertEqual( str(v1) == '[d,u]', True)
        self.assertEqual( v1.describe().replace(' ','') == "-->squark1+[d,u]", True)
        self.assertEqual(v1 > v2, True)  #Larger by number of outgoing particles
        self.assertEqual(v2 > v3, True)  #Larger by mass of sq1
        self.assertEqual(v5b == v5, True)  #All common properties are equal
        
    def testVertexStr(self):
        
        ep = Particle(_name='e+',mass = 0.0*GeV, zParity = 1)
        L = Particle(_name='L', zParity = 1)
        mum = Particle(_name='mu-', eCharge = -1, mass = 0.106*GeV, zParity = 1)
        inP = Particle(zParity = -1)
        outP = Particle(zParity = -1)
        vstr = createVertexFromStr('[e+,L,mu-]',particlesDict)
        v = Vertex(inParticle=inP, outParticles=[ep,mum,L,outP])
        vB = Vertex(inParticle=inP, outParticles=[particlesDict['e+'],
                                                  particlesDict['mu-'],
                                                  particlesDict['ta-'],outP])
        vC = v.copy()
        self.assertEqual(v == vstr, True)  #Check that vertices match
        self.assertEqual(vB == vstr, True)  #Check that vertices match
        self.assertEqual(vB == vC, True)  #Check that vertices match
        for p in vC.outEven:
            if p._name != 'mu-': continue
            self.assertEqual(hasattr(p,'mass'), True)  #Mass has been kept
        
        
        
if __name__ == "__main__":
    unittest.main()
