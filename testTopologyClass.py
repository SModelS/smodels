#!/usr/bin/env python

"""
.. module:: testTopologyClass
   :synopsis: Tests the theory.topology.Topology and TopologyList classes

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.vertex import Vertex
from smodels.theory.branch import Branch
from smodels.theory.element import Element
from smodels.theory.topology import Topology,TopologyList
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
from smodels.tools.physicsUnits import GeV, TeV, fb
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

gluino.mass = 400.*GeV
st1.mass = 300.*GeV
n3.mass = 150.*GeV        
n2.mass = 100.*GeV
n1.mass = 100.*GeV

v0 = Vertex(inParticle=None, outParticles=[gluino])
v0b = Vertex(inParticle=None, outParticles=[st1])
v1 = Vertex(inParticle=gluino, outParticles=[st1,t])
v2 = Vertex(inParticle=st1, outParticles=[n1,bbar,t])
v2b = Vertex(inParticle=st1, outParticles=[n2,bbar,t])
v2c = Vertex(inParticle=st1, outParticles=[n3,bbar,t])

b1 = Branch(vertices=[v0,v1,v2])
b2 = Branch(vertices=[v0b,v2])

b1B = Branch(vertices=[v0,v1,v2b])
b1C = Branch(vertices=[v0,v1,v2c])

el1 = Element(branches=[b1.copy(),b2.copy()])
el1.weight = w1
el1B = Element(branches=[b1B.copy(),b2.copy()])
el1B.weight = w2
el1C = Element(branches=[b1C.copy(),b2.copy()])
el1C.weight = w3     
el2 = Element(branches=[b2.copy(),b2.copy()])
el2.weight = w3

class TopologyTest(unittest.TestCase):
        
    def testTopology(self):
        
        top1 = Topology(elements=[el1.copy()])
        top2 = Topology(elements=[el2.copy()])

        self.assertEqual(str(top1) == "[1,3][1,2,3]", True)
        self.assertEqual(top1 > top2, True)
        self.assertEqual(top2.getElements() == [el2], True)
        self.assertEqual(top1.addElement(el2), False) #Element does not match topology
        self.assertEqual(top2._getTinfo() == {'vertnumb' : [2,2], 'vertparts' : [[1,3],[1,3]]},True)

        top1.addElement(el1B.copy())
        self.assertEqual(len(top1.getElements()) == 1, True)
        self.assertEqual(top1.getElements()[0].weight[0].value == 32.*fb, True)
        self.assertEqual(top1.getTotalWeight()[0].value == 32.*fb, True)
        pids = [[1000006, 1000022], [1000021, 1000006, [1000022, 1000023]]]
        self.assertEqual(top1.getElements()[0].getOddPIDs() == pids, True)
        
        top1.addElement(el1C.copy())
        self.assertEqual(len(top1.getElements()) == 2, True)
        self.assertEqual(top1.getTotalWeight()[0].value == 34.*fb, True)
        self.assertEqual(top1.getElements()[0].weight[0].value == 32.*fb, True)
        self.assertEqual(top1.getElements()[1].weight[0].value == 2.*fb, True)
        pids = [[1000006, 1000022], [1000021, 1000006, 1000025]]
        self.assertEqual(top1.getElements()[1].getOddPIDs() == pids, True)

    def testTopologyList(self):
        
        top1 = Topology(elements=[el1.copy()])
        top2 = Topology(elements=[el2.copy()])
        topL1 = TopologyList(topologies=[top1])
        topL2 = TopologyList(topologies=[top2])
        self.assertEqual(len(topL1) == len(topL2) == 1, True)
        self.assertEqual(topL1 > topL2, True) #Bigger by number of vertices
        
        topL = TopologyList()
        topL.addList(topL1)
        topL.add(top2)
        self.assertEqual(len(topL) == 2, True)
        self.assertEqual(topL.describe() == "TopologyList:\n[1,3][1,3]\n[1,3][1,2,3]\n", True)        
        topL.addElement(el1B.copy())
        topL.addElement(el1C.copy())
        self.assertEqual(topL.getTotalWeight()[0].value == 36.*fb, True)
        
if __name__ == "__main__":
    unittest.main()
