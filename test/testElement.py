#!/usr/bin/env python

"""
.. module:: testElementClass
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.branch import Branch
from smodels.theory.element import Element
from smodels.tools.physicsUnits import GeV, TeV, fb
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
from smodels.share.models import SMparticles, MSSMparticles

from smodels.experiment import finalStateParticles


u = SMparticles.u
d = SMparticles.d
t = SMparticles.t
b = SMparticles.b
g = SMparticles.g
em = SMparticles.e
nue = SMparticles.nue
L = finalStateParticles.LList
e = finalStateParticles.eList


gluino = MSSMparticles.gluino
st1 = MSSMparticles.st1
n1 = MSSMparticles.n1
n2 = MSSMparticles.n2
n3 = MSSMparticles.n3
n4 = MSSMparticles.n4

class ElementTest(unittest.TestCase):
        
    def testElement(self):
         
        b1 = Branch()
        b1.particles = [[t],[b,t]]
        b1.BSMparticles = [gluino,st1,n1]
        b2 = Branch()
        b2.particles = [[b,t]]
        b2.BSMparticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()
         
        el1 = Element()
        el1.branches=[b1,b2]        
        el2 = Element()
        el2.branches=[b2,b2]
        el1B = Element()
        el1B.branches = [b2,b1]
        self.assertEqual(el1 > el2,True) #Bigger by number of vertices
        self.assertFalse(el1 == el1B) #Just differ by branch ordering
        el1.sortBranches()
        el1B.sortBranches()
        self.assertTrue(el1 == el1B) #Now elements should be equal
         
        e1Info = {"vertnumb" : [1,2], "vertparts" : [[2],[1,2]]}
        self.assertEqual(el1.getEinfo() == e1Info, True) 
         
         
    def testElementInclusive(self):
          
        b1 = Branch()
        b1.particles = [[em],[em,nue]]
        b1.BSMparticles = [gluino,st1,n1]
        b1.setInfo()
                 
        b1b = Branch()
        b1b.particles = [[L],[L,nue]]
        b1b.BSMparticles = [gluino,st1,n1]
        b1b.setInfo()
         
        b2 = Branch()
        b2.particles = [[L,nue]]
        b2.BSMparticles = [st1,n1]
        b2.setInfo()
         
        b2b = Branch()
        b2b.particles = [[em,nue]]
        b2b.BSMparticles = [st1,n1]
        b2b.setInfo()
          
        el1 = Element()
        el1.branches = [b1,b2]
        el2 = Element()
        el2.branches = [b2b,b1b]
          
        self.assertFalse(el1 == el2) #Direct comparison should fail
        self.assertTrue(el1.particlesMatch(el2))
         
    def testElementStr(self):
         
        b1 = Branch()
        b1.particles = [[t],[b,t]]
        b1.BSMparticles = [gluino,st1,n1]
        b2 = Branch()
        b2.particles = [[b,t]]
        b2.BSMparticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()
        
        el1 = Element()
        el1.branches=[b1,b2]        
        elstrA = Element('[[[t+],[b,t+]],[[b,t+]]]',finalState=['MET','MET'])
        elstrB = Element('[[[b,t+]],[[t+],[b,t+]]]',finalState=['MET','MET'])
        elstrC = Element('[[[b,t+]],[[t],[b,t]]]',finalState=['MET','MET'])
                 
        self.assertTrue(el1 == elstrA) #Elements should be equal
        self.assertFalse(el1 == elstrB) #Elements should be equal (just switch branches)
        elstrB.sortBranches()
        el1.sortBranches()
        elstrC.sortBranches()
        self.assertTrue(el1 == elstrB) #Elements should be equal (just switch branches)
        self.assertFalse(el1 == elstrC) #Elements should not be identical (distinct labels)
        self.assertTrue(el1.particlesMatch(elstrC)) #Final states should be equal (inclusive labels)
         
    def testElementMassComp(self):
        
        b1 = Branch()
        b1.particles = [[t],[b,t]]
        b1.BSMparticles = [gluino,st1,n1]
        b2 = Branch()
        b2.particles = [[b,t]]
        b2.BSMparticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()
        
        el1 = Element()
        el1.branches=[b1,b2]        
        
         
        #Compress gluino-stop1
        gluino.mass = 400.*GeV
        st1.mass = 398.*GeV
        n1.mass = 390.*GeV        
        el1Comp = el1.massCompress(minmassgap = 5.*GeV)
        b1Comp = Branch()
        b1Comp.particles = [[b,t]]
        b1Comp.BSMparticles = [st1,n1]
        b2Comp = Branch()
        b2Comp.particles = [[b,t]]
        b2Comp.BSMparticles = [st1,n1]
        el2 = Element()
        el2.branches = [b1Comp,b2Comp]
        el2.setEinfo()
        self.assertEqual(el1Comp,el2) #Elements should be equal
        
        
        #Compress stop1-neutralino1
        gluino.mass = 400.*GeV
        st1.mass = 393.*GeV
        n1.mass = 390.*GeV        
        el1Comp = el1.massCompress(minmassgap = 5.*GeV)
  
        b1Comp = Branch()
        b1Comp.particles = [[t]]
        b1Comp.BSMparticles = [gluino,n1]
        b2Comp = Branch()
        b2Comp.particles = []
        b2Comp.BSMparticles = [n1]
        el2 = Element(info=[b1Comp,b2Comp])
        el1.sortBranches()
        el2.sortBranches()
        self.assertEqual(el1Comp,el2) #Elements should be equal        
  
          
        #Compress everything
        el1Comp = el1.massCompress(minmassgap = 10.*GeV) #Fully compress
        b1Comp = Branch()
        b1Comp.particles = []
        b1Comp.BSMparticles = [n1]
        b2Comp = b1Comp.copy()
        el2 = Element(info=[b1Comp,b2Comp])        
        self.assertEqual(el1Comp,el2) #Elements should be equal  
          
 
    def testElementInvComp(self):
           
        gluino.mass = 500.*GeV
        st1.mass = 400.*GeV
        n1.mass = 300.*GeV
        n2.mass = 310.*GeV
        n3.mass = 320.*GeV 
        n4.mass = 330.*GeV  
           
        #Compress one step:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N1 + [nue,nue]
        b1 = Branch()
        b1.particles = [[t],[t],[nue,nue]]
        b1.BSMparticles = [gluino,st1,n3,n1]
        b2 = Branch()
        b2.particles = [[nue,nue]]
        b2.BSMparticles = [n3,n1]
        el1 = Element(info=[b1,b2])
        el1Comp = el1.invisibleCompress()
        
                         
        b1Comp = Branch()
        b1Comp.particles = [[t],[t]]
        b1Comp.BSMparticles = [gluino,st1,n3]
        b2Comp = Branch()
        b2Comp.particles = []
        b2Comp.BSMparticles = [n3]
        el2 = Element(info=[b1Comp,b2Comp])
        el2.sortBranches()
        self.assertEqual(el1Comp,el2) #Elements should be equal
 
        #Compress two steps:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [nue,nue,nue,nue]
        b1 = Branch()
        b1.particles = [[nue,nue]]
        b1.BSMparticles = [n3,n1]        
        b2 = Branch()
        b2.particles = [[t],[t],[nue],[nue,nue,nue,nue]]
        b2.BSMparticles = [gluino,st1,n3,n2,n1]
        el1 = Element(info=[b1,b2])
        
        el1Comp = el1.invisibleCompress()                 
        b1Comp = Branch()
        b1Comp.particles = []
        b1Comp.BSMparticles = [n3]        
        b2Comp = Branch()
        b2Comp.particles = [[t],[t]]
        b2Comp.BSMparticles = [gluino,st1,n3]
        el2 = Element(info=[b1Comp,b2Comp])
        self.assertEqual(el1Comp,el2) #Elements should be equal
         
        #Make sure compression only happens at the end:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [e-,nue,nue,nue]
        b1 = Branch()
        b1.particles = [[nue,nue]]
        b1.BSMparticles = [n3,n1]        
        b2 = Branch()
        b2.particles = [[t],[t],[nue],[e,nue,nue,nue]]
        b2.BSMparticles = [gluino,st1,n3,n2,n1]
        el1 = Element(info=[b1,b2])
        
        el1Comp = el1.invisibleCompress()                 
        b1Comp = Branch()
        b1Comp.particles = []
        b1Comp.BSMparticles = [n3]        
        b2Comp = b2.copy()
        el2 = Element(info=[b1Comp,b2Comp])
        self.assertEqual(el1Comp,el2) #Elements should be equal
 
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
         
          
        b1 = Branch()
        b2 = Branch()
        b1.particles = [[g]]
        b1.BSMparticles = [gluino,n1]
        b2.particles = [[g]]
        b2.BSMparticles = [gluino,n2]
        
        el1 = Element(info=[b1,b1])
        el1.weight = w1
        el2 = Element(info=[b2,b2])
        el2.weight = w2
        el3 = Element(info=[b1,b2])
        el3.weight = w3
        el1.combineWith(el2)
        self.assertEqual(el1.weight[0].value,32.*fb)
        self.assertEqual(el1.getPIDs(),[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        el1.combineWith(el3)
        self.assertEqual(el1.weight[0].value,34.*fb)
        self.assertEqual(el1.getPIDs(),[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        
        
if __name__ == "__main__":
    unittest.main()
