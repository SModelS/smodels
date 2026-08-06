#!/usr/bin/env python3

"""
.. module:: testInclusiveObjects
   :synopsis: Tests implementation of inclusive objects for entries in the database.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.txnameObj import TxName
from smodels.base.physicsUnits import GeV
from smodels.experiment import infoObj
from smodels.base.inclusiveObjects import InclusiveValue, InclusiveList
from smodels.experiment.defaultFinalStates import finalStates
from smodels.experiment.expSMS import ExpSMS
from unitTestHelpers import theorySMSFromString as fromString
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model


slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])

class InclusiveObjectsTest(unittest.TestCase):

    def testInclusiveValue(self):

        x,y = 5,10.
        z = InclusiveValue()
        self.assertTrue(x == z)
        self.assertFalse(y == z)
        self.assertTrue(str(z) == '*')

    def testInclusiveList(self):

        x = [1,'a',10.,-5]
        y = 'string'
        z = InclusiveList()

        self.assertTrue(x == z)
        self.assertFalse(y == z)
        self.assertTrue(str(z) == '[*]')

    def testInclusiveNode(self):

        stringEl = "(PV > gluino(1)), (gluino(1) > st_1(2),e+), (st_1(2) >e-,mu+,N1)"
        sms1 = fromString(stringEl, model=model)
        expInc = ExpSMS.from_string('(PV  > Inclusive(1)), (Inclusive(1) > *anySM,anyBSM)',model=finalStates)

        self.assertTrue(expInc.matchesTo(sms1) is not None)

    def testInclusiveVertex(self):

        stringEl = "[[['q','q'],['e-','nu']],[]]"
        sms1 = fromString(stringEl, model=model, finalState = ['N1','N1'],
                          intermediateState=[['gluino','C1+'],[]])
        expInc = ExpSMS.from_string('(PV  > MET,RHadronG(1)), (RHadronG(1) > *anySM,HSCP+(2)), (HSCP+(2) > *anySM,MET)',
                                    model=finalStates)

        self.assertTrue(expInc.matchesTo(sms1) is not None)

        expIncB = ExpSMS.from_string('(PV  > MET,RHadronG(1)), (RHadronG(1) > *anySM,MET)',
                                    model=finalStates)

        self.assertTrue(expIncB.matchesTo(sms1) is None)


        expIncC = ExpSMS.from_string('(PV  > MET,HSCP+(1)), (HSCP+(1) > RHadronG(2),*anySM), (RHadronG(2) > *anySM,MET)',
                                    model=finalStates)

        self.assertTrue(expIncC.matchesTo(sms1) is None)


    def testInclusiveSMS(self):

        el1 = ExpSMS.from_string("[[['e+'],['e-','mu+']],[[*],['e-','mu+']]]",finalState = ["MET","HSCP"], model=finalStates)
        el2 = ExpSMS.from_string("[[['e+'],['e-','mu+']],[['jet'],['e-','mu+']]]",finalState = ["anyBSM","HSCP"], model=finalStates)
        el3 = ExpSMS.from_string("[[['e+'],['e-','mu+']],[['jet'],['e-','mu+']]]",finalState = ["HSCP","HSCP"], model=finalStates)
        el4 = ExpSMS.from_string("[[*],[['jet'],['e-','mu+']]]",finalState = ["MET","HSCP"], model=finalStates)
        el5 = ExpSMS.from_string("[[['e+'],['e-','mu+']],[['jet'],['e-','mu+']]]",finalState = ["MET","HSCP"], model=finalStates)
        el6 = ExpSMS.from_string("[[*],[*]]",finalState = ["MET","HSCP"], model=finalStates)
        el7 = ExpSMS.from_string("[[['e+'],['e-','mu+']],[['jet','jet'],['e-','mu+']]]",finalState = ["MET","HSCP"], model=finalStates)


        self.assertTrue(el1 == el2)
        self.assertTrue(el2 == el3)
        self.assertTrue(el3 != el4)
        self.assertTrue(el2 == el4)
        self.assertTrue(el1 == el4)
        self.assertTrue(el1 != el3)
        self.assertTrue(el6 == el1)
        self.assertTrue(el6 == el2)
        self.assertTrue(el5 == el1)
        self.assertTrue(el5 != el3)
        self.assertTrue(el7 != el1)
        self.assertTrue(el1 != el7)

        el10 = fromString("[[['q','q'],['e-','nu']],[['e-','nu']]]",
                                finalState = ['N1','N1'],
                                intermediateState=[['gluino','C1+'],['C1+']],
                                model=model)

        el11 = ExpSMS.from_string("(PV > HSCP+(1),HSCP+(2)), (HSCP+(1) > *anySM,MET), (HSCP+(2) > MET,*anySM)",
                                  model=finalStates)

        self.assertFalse(el11 == el10)

        el12 = ExpSMS.from_string("(PV > RHadronG(1),HSCP+(2)), (RHadronG(1) > *anySM,HSCP+(3)), (HSCP+(2) > MET,*anySM), (HSCP+(3) > MET,*anySM)",
                                  model=finalStates)
        self.assertTrue(el12.matchesTo(el10) is not None)

        el13 = fromString("[[['q','q'],['e-','nu']],[]]",
                                finalState = ['N1','N1'], intermediateState=[['gluino','C1+'],[]],
                                model=model)
        el14 = ExpSMS.from_string("(PV > RHadronG(1),HSCP+(2)), (RHadronG(1) > *anySM,MET), (HSCP+(2) > MET,*anySM)",
                                  model=finalStates)

        self.assertFalse(el14 == el13)

        el15 = ExpSMS.from_string("(PV > RHadronG(1),MET), (RHadronG(1) > *anySM,HSCP+(2)), (HSCP+(2) > MET,*anySM)",
                                  model=finalStates)
        self.assertTrue(el15 == el13)

    def testInclusiveTxNameData(self):

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)
        res = tx.txnameData.getValueFor([100.])
        self.assertAlmostEqual(res,0.058038)

        res = tx.txnameData.getValueFor([125.])
        self.assertAlmostEqual(res,0.090999)
        res = tx.txnameData.getValueFor([200])
        self.assertEqual(res,None)

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM6.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)
        res = tx.txnameData.getValueFor([279.,170.,100.])
        self.assertAlmostEqual(res,0.097172,6)
        res = tx.txnameData.getValueFor([1.917E+03,1.7E+02,1E+02])
        self.assertAlmostEqual(res,0.00025745,6)

        res = tx.txnameData.getValueFor([1.112E+03,188,1E+02])
        self.assertAlmostEqual(res,0.015,3)

    def testInclusiveTxName(self):

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)
        sms = fromString("[[],[['e+']]]", intermediateState = [[],['C1+']], finalState = ['C1+','N1'], model=model)

        smsMatch = tx.hasSMSas(sms)  #newEl should be equal to el, but with opposite branch ordering
        self.assertFalse(smsMatch is None)
        nodeDict = dict(zip(smsMatch.nodeIndices,[str(n) for n in smsMatch.nodes]))
        self.assertEqual(nodeDict,{0 : 'PV', 1 : 'C1+', 2 : 'C1+', 3 : 'N1', 4 : 'e+'})


if __name__ == "__main__":
    unittest.main()
