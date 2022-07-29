#!/usr/bin/env python3

"""
.. module:: testElementClass
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from collections import OrderedDict
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
from smodels.share.models import SMparticles
from smodels.experiment.defaultFinalStates import finalStates



slhafile = '../inputFiles/slha/lightEWinos.slha'
model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,erasePrompt=['spin'])

class SMSTest(unittest.TestCase):

    def testTheorySMS(self):


        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(2) >b,t+,N1), (st_1(3) > b,t+,N1)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        sms1 = TheorySMS()
        sms1.add_nodes_from(expSMS.nodes)
        sms1.add_edges_from(expSMS.edgeIndices)
        sms1.prodXSec = 1.0*fb
        sms1.maxWeight = 1.0*fb
        sms1.setGlobalProperties()

        stringEl = "(PV > st_1(1),st_1(2)), (st_1(2) >b,t+,N1), (st_1(1) > b,t+,N1)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        sms2 = TheorySMS()
        sms2.add_nodes_from(expSMS.nodes)
        sms2.add_edges_from(expSMS.edgeIndices)
        sms2.prodXSec = 1.0*fb
        sms2.maxWeight = 1.0*fb
        sms2.setGlobalProperties()

        stringEl = "(PV > st_1(1),gluino(2)), (st_1(1) >b,t+,N1), (gluino(2) > st_1(3),t+), (st_1(3) > b,t+,N1)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        sms1B = TheorySMS()
        sms1B.add_nodes_from(expSMS.nodes)
        sms1B.add_edges_from(expSMS.edgeIndices)
        sms1B.prodXSec = 1.0*fb
        sms1B.maxWeight = 1.0*fb
        sms1B.setGlobalProperties()

        self.assertEqual(sms1 > sms2, True) #Bigger by number of vertices
        self.assertEqual(sms1,sms1B) #Just differ by branch ordering
        self.assertEqual(sms1.canonName, 1110101001101101010000)


    def testExpSMS(self):

        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),e-), (st_1(2) >e-,nue,N1), (st_1(3) > e-,nue,N1)"
        expSMS = ExpSMS.from_string(stringEl, model=model)
        sms1 = TheorySMS()
        sms1.add_nodes_from(expSMS.nodes)
        sms1.add_edges_from(expSMS.edgeIndices)
        sms1.prodXSec = 1.0*fb
        sms1.maxWeight = 1.0*fb
        sms1.setGlobalProperties()

        stringEl = "(PV > RHadronU(1),RHadronG(2)), (RHadronU(1) >e-,nue,MET), (RHadronG(2) > RHadronU(3),L), (RHadronU(3) > L,nue,MET)"
        expSMS2 = ExpSMS.from_string(stringEl, model=finalStates)

        self.assertTrue(expSMS2 == sms1) #SMS match (using inclusive labels)


    def testSMSFromStr(self):

        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(2) >b,t+,N1), (st_1(3) > t+,b,N1)"
        expSMS1 = ExpSMS.from_string(stringEl, model=model)

        expSMSA = ExpSMS.from_string("[[['t+'],['b','t+']],[['b','t+']]]",finalState=['MET','MET'], model=finalStates)
        expSMSB = ExpSMS.from_string("[[['b','t+']],[['t+'],['b','t+']]]",finalState=['MET','MET'], model=finalStates)
        expSMSC = ExpSMS.from_string("[[['b','t+']],[['t'],['b','t']]]",finalState=['MET','MET'], model=finalStates)

        self.assertTrue(expSMS1 == expSMSA) #SMS should be equal
        self.assertTrue(expSMS1 == expSMSB) #SMS should be equal (only branch order differs)
        self.assertTrue(expSMS1 == expSMSC) #SMS should be equal (inclusive labels)

    # def testElementMassComp(self):

    #     b1 = Branch()
    #     b1.evenParticles = [[t],[b,t]]
    #     b1.oddParticles = [gluino,st1,n1]
    #     b2 = Branch()
    #     b2.evenParticles = [[b,t]]
    #     b2.oddParticles = [st1,n1]
    #     b1.setInfo()
    #     b2.setInfo()

    #     el1 = Element()
    #     el1.branches=[b1,b2]


    #     #Compress gluino-stop1
    #     gluino.mass = 400.*GeV
    #     gluino.totalwidth = float('inf')*GeV
    #     st1.mass = 398.*GeV
    #     st1.totalwidth = float('inf')*GeV
    #     n1.mass = 390.*GeV
    #     n1.totalwidth = 0.*GeV
    #     el1Comp = el1.massCompress(minmassgap = 5.*GeV)
    #     b1Comp = Branch()
    #     b1Comp.evenParticles = [[b,t]]
    #     b1Comp.oddParticles = [st1,n1]
    #     b2Comp = Branch()
    #     b2Comp.evenParticles = [[b,t]]
    #     b2Comp.oddParticles = [st1,n1]
    #     el2 = Element()
    #     el2.branches = [b1Comp,b2Comp]
    #     el2.setEinfo()
    #     self.assertEqual(el1Comp,el2) #Elements should be equal


    #     #Compress stop1-neutralino1
    #     gluino.mass = 400.*GeV
    #     st1.mass = 393.*GeV
    #     n1.mass = 390.*GeV
    #     el1Comp = el1.massCompress(minmassgap = 5.*GeV)

    #     b1Comp = Branch()
    #     b1Comp.evenParticles = [[t]]
    #     b1Comp.oddParticles = [gluino,n1]
    #     b2Comp = Branch()
    #     b2Comp.evenParticles = []
    #     b2Comp.oddParticles = [n1]
    #     el2 = Element(info=[b1Comp,b2Comp])
    #     el1.sortBranches()
    #     el2.sortBranches()
    #     self.assertEqual(el1Comp,el2) #Elements should be equal


    #     #Compress everything
    #     el1Comp = el1.massCompress(minmassgap = 10.*GeV) #Fully compress
    #     b1Comp = Branch()
    #     b1Comp.evenParticles = []
    #     b1Comp.oddParticles = [n1]
    #     b2Comp = b1Comp.copy()
    #     el2 = Element(info=[b1Comp,b2Comp])
    #     self.assertEqual(el1Comp,el2) #Elements should be equal


    # def testElementInvComp(self):

    #     gluino.mass = 500.*GeV
    #     st1.mass = 400.*GeV
    #     n1.mass = 300.*GeV
    #     n2.mass = 310.*GeV
    #     n3.mass = 320.*GeV
    #     n4.mass = 330.*GeV

    #     #Compress one step:
    #     #N3 --> N1 + [nue,nue]
    #     #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N1 + [nue,nue]
    #     b1 = Branch()
    #     b1.evenParticles = [[t],[t],[nue,nue]]
    #     b1.oddParticles = [gluino,st1,n3,n1]
    #     b2 = Branch()
    #     b2.evenParticles = [[nue,nue]]
    #     b2.oddParticles = [n3,n1]
    #     el1 = Element(info=[b1,b2])
    #     el1Comp = el1.invisibleCompress()


    #     b1Comp = Branch()
    #     b1Comp.evenParticles = [[t],[t]]
    #     b1Comp.oddParticles = [gluino,st1,n3]
    #     b2Comp = Branch()
    #     b2Comp.evenParticles = []
    #     b2Comp.oddParticles = [n3]
    #     el2 = Element(info=[b1Comp,b2Comp])
    #     el2.sortBranches()
    #     self.assertEqual(el1Comp,el2) #Elements should be equal

    #     #Compress two steps:
    #     #N3 --> N1 + [nue,nue]
    #     #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [nue,nue,nue,nue]
    #     b1 = Branch()
    #     b1.evenParticles = [[nue,nue]]
    #     b1.oddParticles = [n3,n1]
    #     b2 = Branch()
    #     b2.evenParticles = [[t],[t],[nue],[nue,nue,nue,nue]]
    #     b2.oddParticles = [gluino,st1,n3,n2,n1]
    #     el1 = Element(info=[b1,b2])

    #     el1Comp = el1.invisibleCompress()
    #     b1Comp = Branch()
    #     b1Comp.evenParticles = []
    #     b1Comp.oddParticles = [n3]
    #     b2Comp = Branch()
    #     b2Comp.evenParticles = [[t],[t]]
    #     b2Comp.oddParticles = [gluino,st1,n3]
    #     el2 = Element(info=[b1Comp,b2Comp])
    #     self.assertEqual(el1Comp,el2) #Elements should be equal

    #     #Make sure compression only happens at the end:
    #     #N3 --> N1 + [nue,nue]
    #     #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [e-,nue,nue,nue]
    #     b1 = Branch()
    #     b1.evenParticles = [[nue,nue]]
    #     b1.oddParticles = [n3,n1]
    #     b2 = Branch()
    #     b2.evenParticles = [[t],[t],[nue],[e,nue,nue,nue]]
    #     b2.oddParticles = [gluino,st1,n3,n2,n1]
    #     el1 = Element(info=[b1,b2])

    #     el1Comp = el1.invisibleCompress()
    #     b1Comp = Branch()
    #     b1Comp.evenParticles = []
    #     b1Comp.oddParticles = [n3]
    #     b2Comp = b2.copy()
    #     el2 = Element(info=[b1Comp,b2Comp])
    #     self.assertEqual(el1Comp,el2) #Elements should be equal

    # def testElementCombine(self):

    #     gluino.mass = 500.*GeV
    #     st1.mass = 400.*GeV
    #     n1.mass = 250.*GeV
    #     n2.mass = 300.*GeV
    #     n3.mass = 320.*GeV
    #     n1.totalwidth = 0.*GeV
    #     n2.totalwidth = 0.*GeV #just for the sake of the example

    #     w1 = XSectionList()
    #     w1.xSections.append(XSection())
    #     w1.xSections[0].info = XSectionInfo()
    #     w1.xSections[0].info.sqrts = 8.*TeV
    #     w1.xSections[0].info.label = '8 TeV'
    #     w1.xSections[0].info.order = 0
    #     w1.xSections[0].value = 10.*fb
    #     w2 = w1.copy()
    #     w2.xSections[0].value = 22.*fb
    #     w3 = w1.copy()
    #     w3.xSections[0].value = 2.*fb


    #     b1 = Branch()
    #     b2 = Branch()
    #     b1.evenParticles = [[g]]
    #     b1.oddParticles = [gluino,n1]
    #     b2.evenParticles = [[g]]
    #     b2.oddParticles = [gluino,n2]

    #     el1 = Element(info=[b1,b1])
    #     el1.weight = w1
    #     el2 = Element(info=[b2,b2])
    #     el2.weight = w2
    #     el3 = Element(info=[b1,b2])
    #     el3.weight = w3
    #     el1 += el2
    #     self.assertEqual(el1.weight[0].value,32.*fb)
    #     self.assertEqual(el1.pdg,[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
    #     self.assertEqual(el1.getAverage('mass'),[[gluino.mass,(n1.mass+n2.mass)/2.]]*2)
    #     el1 += el3
    #     self.assertEqual(el1.weight[0].value,34.*fb)
    #     self.assertEqual(el1.pdg,[[1000021,[1000022,1000023]],[1000021,[1000022,1000023]]])
    #     self.assertEqual(el1.getAverage('mass'),[[gluino.mass,(n1.mass+n2.mass)/2.]]*2)

if __name__ == "__main__":
    unittest.main()
