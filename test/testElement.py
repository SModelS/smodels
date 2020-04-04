#!/usr/bin/env python3

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
from smodels.share.models import SMparticles, mssm
from smodels.theory.particle import ParticleList

from smodels.experiment.defaultFinalStates import finalStates


u = SMparticles.u
d = SMparticles.d
t = SMparticles.t
b = SMparticles.b
g = SMparticles.g
em = SMparticles.e
nue = SMparticles.nue
L = finalStates.getParticlesWith(label='L')[0]
e = finalStates.getParticlesWith(label='e')[0]


gluino = mssm.gluino
st1 = mssm.st1
n1 = mssm.n1
n2 = mssm.n2
n3 = mssm.n3
n4 = mssm.n4

class ElementTest(unittest.TestCase):

    def testElement(self):

        b1 = Branch()
        b1.evenParticles = [[t],[b,t]]
        b1.oddParticles = [gluino,st1,n1]
        b2 = Branch()
        b2.evenParticles = [[b,t]]
        b2.oddParticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()

        el1 = Element()
        el1.branches=[b1,b2]
        el2 = Element()
        el2.branches=[b2,b2]
        el1B = Element()
        el1B.branches = [b2,b1]
        self.assertEqual(el1 > el2,True) #Bigger by number of vertices
        self.assertEqual(el1,el1B) #Just differ by branch ordering
        el1.sortBranches()
        e1Info = {"vertnumb" : [1,2], "vertparts" : [[2],[1,2]]}
        self.assertEqual(el1.getEinfo() == e1Info, True)


    def testElementInclusive(self):

        b1 = Branch()
        b1.evenParticles = [[em],[em,nue]]
        b1.oddParticles = [gluino,st1,n1]
        b1.setInfo()

        b1b = Branch()
        b1b.evenParticles = [[L],[L,nue]]
        b1b.oddParticles = [gluino,st1,n1]
        b1b.setInfo()

        b2 = Branch()
        b2.evenParticles = [[L,nue]]
        b2.oddParticles = [st1,n1]
        b2.setInfo()

        b2b = Branch()
        b2b.evenParticles = [[em,nue]]
        b2b.oddParticles = [st1,n1]
        b2b.setInfo()

        el1 = Element()
        el1.branches = [b1,b2]
        el2 = Element()
        el2.branches = [b2b,b1b]

        self.assertTrue(el1 == el2) #Elements match (using inclusive labels)

    def testElementStr(self):

        b1 = Branch()
        b1.evenParticles = [ParticleList([t]),ParticleList([b,t])]
        b1.oddParticles = [gluino,st1,n1]
        b2 = Branch()
        b2.evenParticles = [ParticleList([b,t])]
        b2.oddParticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()

        el1 = Element()
        el1.branches=[b1,b2]
        elstrA = Element('[[[t+],[b,t+]],[[b,t+]]]',finalState=['MET','MET'], model=finalStates)
        elstrB = Element('[[[b,t+]],[[t+],[b,t+]]]',finalState=['MET','MET'], model=finalStates)
        elstrC = Element('[[[b,t+]],[[t],[b,t]]]',finalState=['MET','MET'], model=finalStates)


        self.assertTrue(el1 == elstrA) #Elements should be equal
        self.assertTrue(el1 == elstrB) #Elements should be equal (only branch order differs)
        self.assertTrue(el1 == elstrC) #Elements should be equal (inclusive labels)

    def testElementMassComp(self):

        b1 = Branch()
        b1.evenParticles = [[t],[b,t]]
        b1.oddParticles = [gluino,st1,n1]
        b2 = Branch()
        b2.evenParticles = [[b,t]]
        b2.oddParticles = [st1,n1]
        b1.setInfo()
        b2.setInfo()

        el1 = Element()
        el1.branches=[b1,b2]


        #Compress gluino-stop1
        gluino.mass = 400.*GeV
        gluino.totalwidth = float('inf')*GeV
        st1.mass = 398.*GeV
        st1.totalwidth = float('inf')*GeV
        n1.mass = 390.*GeV
        n1.totalwidth = 0.*GeV
        el1Comp = el1.massCompress(minmassgap = 5.*GeV)
        b1Comp = Branch()
        b1Comp.evenParticles = [[b,t]]
        b1Comp.oddParticles = [st1,n1]
        b2Comp = Branch()
        b2Comp.evenParticles = [[b,t]]
        b2Comp.oddParticles = [st1,n1]
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
        b1Comp.evenParticles = [[t]]
        b1Comp.oddParticles = [gluino,n1]
        b2Comp = Branch()
        b2Comp.evenParticles = []
        b2Comp.oddParticles = [n1]
        el2 = Element(info=[b1Comp,b2Comp])
        el1.sortBranches()
        el2.sortBranches()
        self.assertEqual(el1Comp,el2) #Elements should be equal


        #Compress everything
        el1Comp = el1.massCompress(minmassgap = 10.*GeV) #Fully compress
        b1Comp = Branch()
        b1Comp.evenParticles = []
        b1Comp.oddParticles = [n1]
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
        b1.evenParticles = [[t],[t],[nue,nue]]
        b1.oddParticles = [gluino,st1,n3,n1]
        b2 = Branch()
        b2.evenParticles = [[nue,nue]]
        b2.oddParticles = [n3,n1]
        el1 = Element(info=[b1,b2])
        el1Comp = el1.invisibleCompress()


        b1Comp = Branch()
        b1Comp.evenParticles = [[t],[t]]
        b1Comp.oddParticles = [gluino,st1,n3]
        b2Comp = Branch()
        b2Comp.evenParticles = []
        b2Comp.oddParticles = [n3]
        el2 = Element(info=[b1Comp,b2Comp])
        el2.sortBranches()
        self.assertEqual(el1Comp,el2) #Elements should be equal

        #Compress two steps:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [nue,nue,nue,nue]
        b1 = Branch()
        b1.evenParticles = [[nue,nue]]
        b1.oddParticles = [n3,n1]
        b2 = Branch()
        b2.evenParticles = [[t],[t],[nue],[nue,nue,nue,nue]]
        b2.oddParticles = [gluino,st1,n3,n2,n1]
        el1 = Element(info=[b1,b2])

        el1Comp = el1.invisibleCompress()
        b1Comp = Branch()
        b1Comp.evenParticles = []
        b1Comp.oddParticles = [n3]
        b2Comp = Branch()
        b2Comp.evenParticles = [[t],[t]]
        b2Comp.oddParticles = [gluino,st1,n3]
        el2 = Element(info=[b1Comp,b2Comp])
        self.assertEqual(el1Comp,el2) #Elements should be equal

        #Make sure compression only happens at the end:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [e-,nue,nue,nue]
        b1 = Branch()
        b1.evenParticles = [[nue,nue]]
        b1.oddParticles = [n3,n1]
        b2 = Branch()
        b2.evenParticles = [[t],[t],[nue],[e,nue,nue,nue]]
        b2.oddParticles = [gluino,st1,n3,n2,n1]
        el1 = Element(info=[b1,b2])

        el1Comp = el1.invisibleCompress()
        b1Comp = Branch()
        b1Comp.evenParticles = []
        b1Comp.oddParticles = [n3]
        b2Comp = b2.copy()
        el2 = Element(info=[b1Comp,b2Comp])
        self.assertEqual(el1Comp,el2) #Elements should be equal

    def testElementCombine(self):

        gluino.mass = 500.*GeV
        st1.mass = 400.*GeV
        n1.mass = 250.*GeV
        n2.mass = 300.*GeV
        n3.mass = 320.*GeV
        n1.totalwidth = 0.*GeV
        n2.totalwidth = 0.*GeV #just for the sake of the example

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
        b1.evenParticles = [[g]]
        b1.oddParticles = [gluino,n1]
        b2.evenParticles = [[g]]
        b2.oddParticles = [gluino,n2]

        el1 = Element(info=[b1,b1])
        el1.weight = w1
        el2 = Element(info=[b2,b2])
        el2.weight = w2
        el3 = Element(info=[b1,b2])
        el3.weight = w3
        el1 += el2
        self.assertEqual(el1.weight[0].value,32.*fb)
        self.assertEqual(el1.pdg,[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        self.assertEqual(el1.getAverage('mass'),[[gluino.mass,(n1.mass+n2.mass)/2.]]*2)
        el1 += el3
        self.assertEqual(el1.weight[0].value,34.*fb)
        self.assertEqual(el1.pdg,[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
        self.assertEqual(el1.getAverage('mass'),[[gluino.mass,(n1.mass+n2.mass)/2.]]*2)

if __name__ == "__main__":
    unittest.main()
