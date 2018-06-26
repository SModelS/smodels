#!/usr/bin/env python

"""
.. module:: testElementClass
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.vertex import Vertex
from smodels.theory.branch import Branch
from smodels.theory.element import Element,createElementFromStr
from smodels.tools.physicsUnits import GeV, TeV, fb
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
import pickle

#Load the particle dictionaries
f = open("particleDefinitions.pcl","rb")
modelParticles = pickle.load(f)
particlesDict = dict([[p._name,p] for p in modelParticles])
f.close()

u = particlesDict['u']
d = particlesDict['d']
t = particlesDict['t+']
bbar = particlesDict['b+']
g = particlesDict['g']
em = particlesDict['e-']
nue = particlesDict['nue']
L = particlesDict['L']
e = particlesDict['e']


gluino = particlesDict['gluino']
st1 = particlesDict['st_1']
n1 = particlesDict['N1']
n2 = particlesDict['N2']
n3 = particlesDict['N3']
n4 = particlesDict['N4']

class ElementTest(unittest.TestCase):
        
    def testElement(self):
        
        v0 = Vertex(inParticle=None, outParticles=[gluino])
        v0b = Vertex(inParticle=None, outParticles=[st1])
        v1 = Vertex(inParticle=gluino, outParticles=[st1,t])
        v2 = Vertex(inParticle=st1, outParticles=[n1,bbar,t])
        
        b1 = Branch(vertices=[v0,v1,v2])
        b2 = Branch(vertices=[v0b,v2])
        
        el1 = Element(branches=[b1,b2])        
        el2 = Element(branches=[b2,b2])
        el1B = Element()
        el1B.branches = [b2,b1]
        self.assertEqual(el1 > el2,True) #Bigger by number of vertices
        self.assertEqual(el1 == el1B,True) #Just differ br branch ordering
        
        e1Info = {"vertnumb" : [2,3], "vertparts" : [[1,3],[1,2,3]]}
        self.assertEqual(el1.getEinfo() == e1Info, True) 
        
        
    def testElementInclusive(self):
        
        v0 = Vertex(inParticle=None, outParticles=[gluino])
        v0b = Vertex(inParticle=None, outParticles=[st1])     
        v1 = Vertex(inParticle=gluino, outParticles=[st1,em])
        v2 = Vertex(inParticle=st1, outParticles=[n1,em,nue])

        v1b = Vertex(inParticle=gluino, outParticles=[st1,L])
        v2b = Vertex(inParticle=st1, outParticles=[n1,L,nue])
        
        b1 = Branch(vertices = [v0,v1,v2])
        b1b = Branch(vertices = [v0,v1b,v2b])
        b2 = Branch(vertices=[v0b,v2b])
        b2b = Branch(vertices=[v0b,v2])
        
        el1 = Element(branches = [b1,b2])
        el2 = Element()
        el2.branches = [b2b,b1b]
        
        self.assertEqual( el1 == el2, True) #Test if inclusive label comparison works
        
    def testElementStr(self):
        
        v0 = Vertex(inParticle=None, outParticles=[gluino])
        v0b = Vertex(inParticle=None, outParticles=[st1])
        v1 = Vertex(inParticle=gluino, outParticles=[st1,t])
        v2 = Vertex(inParticle=st1, outParticles=[n1,bbar,t])
        
        b1 = Branch(vertices=[v0,v1,v2])
        b2 = Branch(vertices=[v0b,v2])
        
        el1 = Element(branches=[b1,b2])
                        
        elstrA = createElementFromStr('[[[t+],[b+,t+]],[[b+,t+]]]',particlesDict)
        elstrB = createElementFromStr('[[[b+,t+]],[[t+],[b+,t+]]]',particlesDict)
        elstrC = createElementFromStr('[[[b,t+]],[[t],[b+,t]]]',particlesDict)
                
        self.assertEqual( el1 == elstrA, True) #Elements should be equal
        self.assertEqual( el1 == elstrB, True) #Elements should be equal (just switch branches)
        self.assertEqual( el1 == elstrC, True) #Elements should be equal (inclusive labels)
        
    def testElementMassComp(self):


        v0 = Vertex(inParticle=None, outParticles=[gluino])
        v0b = Vertex(inParticle=None, outParticles=[st1])
        v0c = Vertex(inParticle=None, outParticles=[n1])
        v1 = Vertex(inParticle=gluino, outParticles=[st1,t])
        v2 = Vertex(inParticle=st1, outParticles=[n1,bbar,t])
        v1c = Vertex(inParticle=gluino, outParticles=[n1,t])
        
        b1 = Branch(vertices=[v0,v1,v2])
        b2 = Branch(vertices=[v0b,v2])
        el1 = Element(branches=[b1,b2])
        
        #Compress gluino-stop1
        gluino.mass = 400.*GeV
        st1.mass = 398.*GeV
        n1.mass = 390.*GeV        
        el1Comp = el1.massCompress(minmassgap = 5.*GeV)
        b1Comp = Branch(vertices=[v0b,v2])
        b2Comp = b2
        el2 = Element(branches=[b1Comp,b2Comp])
        self.assertEqual( el1Comp == el2, True) #Elements should be equal
        
        
        #Compress stop1-neutralino1
        gluino.mass = 400.*GeV
        st1.mass = 393.*GeV
        n1.mass = 390.*GeV        
        el1Comp = el1.massCompress(minmassgap = 5.*GeV)

        b1Comp = Branch(vertices=[v0,v1c])
        b2Comp = Branch(vertices=[v0c])
        el2 = Element(branches=[b1Comp,b2Comp])
        self.assertEqual( el1Comp == el2, True) #Elements should be equal        

        
        #Compress everything
        el1Comp = el1.massCompress(minmassgap = 10.*GeV) #Fully compress
        b1Comp = Branch(vertices=[v0c])
        b2Comp = Branch(vertices=[v0c])
        el2 = Element(branches=[b1Comp,b2Comp])        
        self.assertEqual( el1Comp == el2, True) #Elements should be equal  
        

    def testElementInvComp(self):
          
        gluino.mass = 500.*GeV
        st1.mass = 400.*GeV
        n1.mass = 300.*GeV
        n2.mass = 310.*GeV
        n3.mass = 320.*GeV 
        n4.mass = 330.*GeV  
          
        v0 = Vertex(inParticle=None, outParticles=[gluino])        
        v0c = Vertex(inParticle=None, outParticles=[n3])
        v1 = Vertex(inParticle=gluino, outParticles=[st1,t])
        v2 = Vertex(inParticle=st1, outParticles=[n3,t])
        v3 = Vertex(inParticle=n3, outParticles=[n1,nue,nue])        
        v4 = Vertex(inParticle=n3, outParticles=[n2,nue])
        v5 = Vertex(inParticle=n2, outParticles=[n1,nue,nue,nue,nue])
        v6 = Vertex(inParticle=n2, outParticles=[n1,nue,nue,nue,em])
          

         
        #Compress one step:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N1 + [nue,nue]
        b1 = Branch(vertices=[v0,v1,v2,v3])
        b2 = Branch(vertices=[v0c,v3])
        el1 = Element(branches=[b1,b2])
        el1Comp = el1.invisibleCompress()                 
        b1Comp = Branch(vertices=[v0,v1,v2])
        b2Comp = Branch(vertices=[v0c])
        el2 = Element(branches=[b1Comp,b2Comp])
        self.assertEqual( el1Comp == el2, True) #Elements should be equal

        #Compress two steps:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [nue,nue,nue,nue]
        b1 = Branch(vertices=[v0,v1,v2,v4,v5])
        b2 = Branch(vertices=[v0c,v3])
        el1 = Element(branches=[b1,b2])
        el1Comp = el1.invisibleCompress()                 
        b1Comp = Branch(vertices=[v0,v1,v2])
        b2Comp = Branch(vertices=[v0c])
        el2 = Element(branches=[b1Comp,b2Comp])
        self.assertEqual( el1Comp == el2, True) #Elements should be equal
        
        #Make sure compression only happens at the end:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [e-,nue,nue,nue]
        b1 = Branch(vertices=[v0,v1,v2,v4,v6])
        b2 = Branch(vertices=[v0c,v3])
        el1 = Element(branches=[b1,b2])
        el1Comp = el1.invisibleCompress()                 
        b1Comp = Branch(vertices=[v0,v1,v2,v4,v6])
        b2Comp = Branch(vertices=[v0c])
        el2 = Element(branches=[b1Comp,b2Comp])
        self.assertEqual( el1Comp == el2, True) #Elements should be equal

    def testElementCombine(self):
        
        gluino.mass = 500.*GeV
        st1.mass = 400.*GeV
        n1.mass = 300.*GeV
        n2.mass = 300.*GeV
        n3.mass = 320.*GeV 
        
        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb
        w2 = w1.copy()
        w2.xSections[0].value = 22.*fb
        w3 = w1.copy()
        w3.xSections[0].value = 2.*fb
        
         
        v0 = Vertex(inParticle=None, outParticles=[gluino])
        v1 = Vertex(inParticle=gluino, outParticles=[g,n1])
        v2 = Vertex(inParticle=gluino, outParticles=[g,n2])
         
        b1 = Branch(vertices=[v0,v1])
        b2 = Branch(vertices=[v0,v2])
        el1 = Element(branches=[b1.copy(),b1.copy()])
        el1.weight = w1
        el2 = Element(branches=[b2.copy(),b2.copy()])
        el2.weight = w2
        el3 = Element(branches=[b1.copy(),b2.copy()])
        el3.weight = w3
        el1.combineWith(el2)
        self.assertEqual(el1.weight[0].value,32.*fb)
        self.assertEqual(el1.getOddPIDs(),[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        el1.combineWith(el3)
        self.assertEqual(el1.weight[0].value,34.*fb)
        self.assertEqual(el1.getOddPIDs(),[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        
                
        
if __name__ == "__main__":
    unittest.main()
