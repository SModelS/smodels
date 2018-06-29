#!/usr/bin/env python

"""
.. module:: testTopologyClass
   :synopsis: Tests the theory.topology.Topology and TopologyList classes

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models import SMparticles, MSSMparticles
from smodels.theory.branch import Branch
from smodels.theory.element import Element
from smodels.theory.topology import Topology,TopologyList
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
from smodels.theory.model import Model
from smodels.tools.physicsUnits import GeV, TeV, fb
import pickle



u = SMparticles.u
d = SMparticles.d
t = SMparticles.t
b = SMparticles.b
g = SMparticles.g
em = SMparticles.e
nue = SMparticles.nue

gluino = MSSMparticles.gluino
st1 = MSSMparticles.st1
n1 = MSSMparticles.n1
n2 = MSSMparticles.n2
n3 = MSSMparticles.n3
n4 = MSSMparticles.n4

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
"""
setattr(getParticlesWith(label = 'gluino'), 'mass', 400.*GeV)
setattr(getParticlesWith(label = 'st1'), 'mass', 300.*GeV)
setattr(getParticlesWith(label = 'n3'), 'mass', 150.*GeV)
setattr(getParticlesWith(label = 'n2'), 'mass', 100.*GeV)
setattr(getParticlesWith(label = 'n1'), 'mass', 100.*GeV)
"""
b1 = Branch()
b1.particles = [[t],[b,t]]
b1.BSMparticles = [gluino,st1,n1]
b1b = Branch()
b1b.particles = [[t],[b,t]]
b1b.BSMparticles = [gluino,st1,n2]
b2 = Branch()
b2.particles = [[b,t]]
b2.BSMparticles = [st1,n1]
b1.setInfo()
b2.setInfo()
 
el1 = Element()
el1.branches=[b1,b2] 
el1.weight = w1
el2 = Element()
el2.branches=[b2,b2]
el2.weight = w3
el1B = Element()
el1B.branches = [b1b,b2]
el1B.weight = w2


class TopologyTest(unittest.TestCase):
        
    def testTopology(self):
        
        top1 = Topology(elements=[el1.copy()])
        top2 = Topology(elements=[el2.copy()])

        self.assertEqual(str(top1) == "[1,2][2]", True)
        self.assertEqual(top1 > top2, True)
        self.assertEqual(top2.getElements() == [el2], True)
        self.assertEqual(top1.addElement(el2), False) #Element does not match topology
        self.assertEqual(top2._getTinfo() == {'vertnumb' : [1,1], 'vertparts' : [[2],[2]]},True)

        top1.addElement(el1B.copy())
        self.assertEqual(len(top1.getElements()) == 1, True)
        self.assertEqual(top1.getElements()[0].weight[0].value == 32.*fb, True)
        self.assertEqual(top1.getTotalWeight()[0].value == 32.*fb, True)
        bsmparticles = [[gluino,st1,n1], [st1,n1]]
        self.assertEqual(top1.getElements()[0].getBSMparticles() == bsmparticles, True)


    def testTopologyList(self):
        
        top1 = Topology(elements=[el1.copy()])
        top2 = Topology(elements=[el2.copy()])
        topL1 = TopologyList(topologies=[top1])
        topL2 = TopologyList(topologies=[top2])
        
        self.assertEqual(len(topL1) == len(topL2) == 1, True)
        self.assertEqual(top1 < top2, False) #Bigger by number of vertices
        
        topL = TopologyList()
        topL.addList(topL1)
        topL.add(top2)
        self.assertEqual(len(topL) == 2, True)
        self.assertEqual(topL.describe() == "TopologyList:\n[2][2]\n[1,2][2]\n", True)        
        topL.addElement(el1B.copy())
        self.assertEqual(topL.getTotalWeight()[0].value == 34.*fb, True)
        
if __name__ == "__main__":
    unittest.main()
