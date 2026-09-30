#!/usr/bin/env python3

"""
.. module:: testElementClass
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.particle import Particle
from smodels.base.model import Model
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.base.crossSection import XSection,XSectionInfo,XSectionList
from smodels.experiment.defaultFinalStates import finalStates
from unitTestHelpers import theorySMSFromString as fromString



slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])
invisible = Particle(label='invisible',pdg=50000,mass=500*GeV,isSM=False)
model.BSMparticles.append(invisible)

class SMSTest(unittest.TestCase):

    def testTheorySMS(self):


        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(2) >b,t+,N1), (st_1(3) > b,t+,N1)"
        sms1 = fromString(stringEl, model=model)

        stringEl = "(PV > st_1(1),st_1(2)), (st_1(2) >b,t+,N1), (st_1(1) > b,t+,N1)"
        sms2 = fromString(stringEl, model=model)

        stringEl = "(PV > st_1(1),gluino(2)), (st_1(1) >b,t+,N1), (gluino(2) > st_1(3),t+), (st_1(3) > b,t+,N1)"
        sms1B = fromString(stringEl, model=model)

        self.assertEqual(sms1 > sms2, True) #Bigger by number of vertices
        self.assertEqual(sms1,sms1B) #Just differ by branch ordering
        self.assertEqual(sms1.canonName, 1110101001101101010000)


    def testExpSMS(self):

        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),e-), (st_1(2) >e-,nue,N1), (st_1(3) > e-,nue,N1)"
        sms1 = fromString(stringEl, model=model)

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

    def testSMSMassComp(self):


        stringEl = "(PV > gluino(1),st_1(2)), (gluino(1) > st_1(3),t+), (st_1(2) >b,t+,N1), (st_1(3) > b,t+,N1)"
        sms1 = fromString(stringEl, model=model)

        gluino = model.getParticle(label='gluino')
        st1 = model.getParticle(label='st_1')
        n1 = model.getParticle(label='N1')

        #Compress gluino-stop1
        gluino.mass = 400.*GeV
        gluino.totalwidth = float('inf')*GeV
        st1.mass = 398.*GeV
        st1.totalwidth = float('inf')*GeV
        n1.mass = 390.*GeV
        n1.totalwidth = 0.*GeV

        smsComp = sms1.massCompress(minmassgap = 5.*GeV, minmassgapISR= 5*GeV)

        stringEl = "(PV > st_1(1),st_1(2)), (st_1(2) >b,t+,N1), (st_1(1) > b,t+,N1)"
        sms2 = fromString(stringEl, model=model)

        self.assertEqual(smsComp,sms2) #Elements should be equal


        #Compress stop1-neutralino1
        gluino.mass = 400.*GeV
        st1.mass = 393.*GeV
        n1.mass = 390.*GeV
        smsComp = sms1.massCompress(minmassgap = 5.*GeV, minmassgapISR= 5*GeV)

        stringEl = "(PV > gluino(1),N1), (gluino(1) > t+,N1)"
        sms2 = fromString(stringEl, model=model)

        self.assertEqual(smsComp,sms2) #Elements should be equal


        #Compress everything
        smsComp = sms1.massCompress(minmassgap = 10.*GeV, minmassgapISR=10*GeV) #Fully compress
        stringEl = "(PV > N1,N1)"
        sms2 = fromString(stringEl, model=model)

        self.assertEqual(smsComp,sms2) #Elements should be equal


    def testSMSInvComp(self):


        gluino = model.getParticle(label='gluino')
        st1 = model.getParticle(label='st_1')
        n1 = model.getParticle(label='N1')
        n2 = model.getParticle(label='N2')
        n3 = model.getParticle(label='N3')
        n4 = model.getParticle(label='N4')
        inv = model.getParticle(label='invisible')


        gluino.mass = 500.*GeV
        st1.mass = 400.*GeV
        n1.mass = 300.*GeV
        n2.mass = 310.*GeV
        n3.mass = 320.*GeV
        n4.mass = 330.*GeV
        inv.mass = n3.mass

        #Compress one step:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N1 + [nue,nue]
        stringEl = "(PV > gluino(1),N3(2)), (gluino(1) > st_1(3),t+), (N3(2) >nue,nue,N1), (st_1(3) > N3(4),t+), (N3(4) > nue,nue,N1)"
        sms1 = fromString(stringEl, model=model)
        smsComp = sms1.invisibleCompress()

        stringEl = "(PV > gluino(1),invisible), (gluino(1) > st_1(3),t+), (st_1(3) > invisible,t+)"
        sms2 = fromString(stringEl, model=model)
        self.assertEqual(smsComp,sms2) #SMS should be equal

        #Compress two steps:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [nue,nue,nue,nue]
        stringEl = "(PV > N3(2),gluino(1)), (N3(2) > nue,nue,N1), (gluino(1) > st_1(3),t+), (st_1(3) > N3(4),t+), (N3(4) > nue,N2(5)), (N2(5)>nue,nue,nue,nue,N1)"
        sms1 = fromString(stringEl, model=model)
        smsComp = sms1.invisibleCompress()

        stringEl = "(PV > invisible,gluino(1)), (gluino(1) > st_1(3),t+), (st_1(3) > invisible,t+)"
        sms2 = fromString(stringEl, model=model)
        self.assertEqual(smsComp,sms2) # SMS should be equal

        #Make sure compression only happens at the end:
        #N3 --> N1 + [nue,nue]
        #gluino --> st_1 + [t+]/st_1 --> N3 + [t+]/N3 --> N2 + [nue]/N2 --> N1 + [e-,nue,nue,nue]
        stringEl = "(PV > N3(2),gluino(1)), (N3(2) > nue,nue,N1), (gluino(1) > st_1(3),t+), (st_1(3) > N3(4),t+), (N3(4) > nue,N2(5)), (N2(5)>e-,nue,nue,nue,N1)"
        sms1 = fromString(stringEl, model=model)
        smsComp = sms1.invisibleCompress()

        stringEl = "(PV > invisible,gluino(1)), (gluino(1) > st_1(3),t+), (st_1(3) > N3(4),t+), (N3(4) > nue,N2(5)), (N2(5)>e-,nue,nue,nue,N1)"
        sms2 = fromString(stringEl, model=model)
        self.assertEqual(smsComp,sms2) #SMS should be equal


    def testElementCombine(self):

        gluino = model.getParticle(label='gluino')
        st1 = model.getParticle(label='st_1')
        n1 = model.getParticle(label='N1')
        n2 = model.getParticle(label='N2')
        n3 = model.getParticle(label='N3')

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


        stringEl = "(PV > gluino(2),gluino(1)), (gluino(2) > g,N1), (gluino(1) > N1,g)"
        sms1 = fromString(stringEl, model=model,prodXSec = w1, maxWeight = w1.getMaxXsec())

        stringEl = "(PV > gluino(2),gluino(1)), (gluino(2) > g,N2), (gluino(1) > N2,g)"
        sms2 = fromString(stringEl, model=model,prodXSec = w2, maxWeight = w2.getMaxXsec())

        stringEl = "(PV > gluino(2),gluino(1)), (gluino(2) > g,N1), (gluino(1) > N2,g)"
        sms3 = fromString(stringEl, model=model,prodXSec = w3, maxWeight = w3.getMaxXsec())

        sms1 += sms2
        self.assertEqual(sms1.weightList[0].value,32.*fb)
        self.assertEqual(sms1.pdg,[None, 1000021, 1000021,
                                    [1000022, 1000023], [21, -21],
                                    [1000022, 1000023], [21, -21]])

        sms1 += sms3
        self.assertEqual(sms1.weightList[0].value,34.*fb)
        self.assertEqual(sms1.pdg,[None, 1000021, 1000021,
                                    [1000022, 1000023], [21, -21],
                                    [1000022, 1000023], [21, -21]])


if __name__ == "__main__":
    unittest.main()
