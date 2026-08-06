#!/usr/bin/env python3

"""
.. module:: testTopologyClass
   :synopsis: Tests the theory.topology.Topology and TopologyList classes

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models import SMparticles, mssm
from smodels.decomposition.topologyDict import TopologyDict
from smodels.base.crossSection import XSection,XSectionInfo,XSectionList
from smodels.base.physicsUnits import TeV, fb

from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.particle import Particle
from smodels.base.model import Model
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.experiment.defaultFinalStates import finalStates
from unitTestHelpers import theorySMSFromString as fromString



slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])
invisible = Particle(label='invisible',pdg=50000,mass=500*GeV,isSM=False)
model.BSMparticles.append(invisible)


u = model.getParticle(label='u')
d = model.getParticle(label='d')
t = model.getParticle(label='t+')
b = model.getParticle(label='b')
g = model.getParticle(label='g')
em = model.getParticle(label='e-')
nue = model.getParticle(label='nue')

gluino = model.getParticle(label='gluino')
st1 = model.getParticle(label='st_1')
n1 = model.getParticle(label='N1')
n2 = model.getParticle(label='N2')
n3 = model.getParticle(label='N3')
n4 = model.getParticle(label='N4')

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

sms1 = fromString('(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(3) > b,t+,N1), (st_1(2) > b,t+,N1)',
                  model=model,
                  prodXSec=w1)
sms2 = fromString('(PV > st_1(1),st_1(2)), (st_1(1) > b,t+,N1), (st_1(2) > b,t+,N1)',
                  model=model,
                  prodXSec=w3)
sms1B = fromString('(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(3) > b,t+,N1), (st_1(2) > b,t+,N1)',
                   model=model,
                   prodXSec=w2)


class TopologyTest(unittest.TestCase):
        
    def testTopologyDict(self):
        
        top1 = TopologyDict()
        top1.addSMS(sms1.copy())
        top2 = TopologyDict()
        top2.addSMS(sms2.copy())

        self.assertEqual(list(top1.keys()),[1110101001101101010000])
        self.assertEqual(list(top2.keys()), [111010100110101000])
        self.assertEqual(top2.getSMSList(),[sms2])
        top1.addSMS(sms2)
        self.assertEqual(len(top1), 2)

        top1.addSMS(sms1B)
        self.assertEqual(len(top1.getSMSList()), 2)
        self.assertEqual(len(top1.getSMSList(1110101001101101010000)), 1)
        self.assertEqual(top1.getSMSList(1110101001101101010000)[0].weightList[0].value,32.*fb)
        self.assertEqual(top1.getTotalWeight()[0].value == 34.*fb, True)
        top1sms1 = top1.getSMSList(1110101001101101010000)[0]
        bsmParticles = [node.particle for node in top1sms1.nodes
                        if node.particle.isSM is False]
        # bsmparticles = [[gluino,st1,n1], [st1,n1]]
        self.assertEqual(bsmParticles,[st1,gluino,n1,st1,n1])

        
if __name__ == "__main__":
    unittest.main()
