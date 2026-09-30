#!/usr/bin/env python3

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from smodels.share.models import mssm
from smodels.experiment import infoObj
from smodels.experiment.txnameObj import TxName, rescaleWidth, unscaleWidth
from smodels.base.physicsUnits import GeV
from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.experiment.expAuxiliaryFuncs import flattenArray
from databaseLoader import database
from smodels.experiment.defaultFinalStates import finalStates
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from unitTestHelpers import theorySMSFromString as fromString
import numpy as np
import unittest

c1 = mssm.c1
n1 = mssm.n1

class TxTest(unittest.TestCase):
    def testWidthTransformation(self):
        a = 1e-10*GeV
        b = unscaleWidth(rescaleWidth(a))*GeV
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 1e-50*GeV
        b = unscaleWidth(rescaleWidth(a))*GeV
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = float('inf')*GeV
        b = unscaleWidth(rescaleWidth(a))*GeV
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 1e-25*GeV
        b = unscaleWidth(rescaleWidth(a))*GeV
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 0.*GeV
        b = unscaleWidth(rescaleWidth(a))*GeV
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))

    def testTxnameSMS(self):

        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        self.assertTrue(len(tx.smsMap), 1)
        sms = list(tx.smsMap.keys())[0]
        defaultSMS = '(PV > Inclusive(1),HSCP), (Inclusive(1) > MET,*anySM)'
        self.assertEqual(str(sms).replace(' ',''), defaultSMS.replace(' ',''))

    def testBrokenFinalState(self):

        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2broken.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        try:
            #print (gInfo)
            TxName(f,gInfo,gInfo,finalStates)
            #print (tx)
            gotError = False
        except SModelSError as e:
            gotError = e

        errstr = "has not been found"
        self.assertTrue(errstr in gotError.args[0])
        self.assertTrue('non-MET' in gotError.args[0])

    def testgetEffFor(self):

        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model( BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])
        sms = fromString("[[],[]]",finalState = ['C1+','N1'],  model=model)

        n1 = model.getParticle(label='N1')
        c1 = model.getParticle(label='C1+')
        n1.mass = 200.0*GeV
        c1.mass = 150.0*GeV
        n1.totalwidth = 0.0*GeV
        c1.totalwidth = 1e-30*GeV

        smsMatch = tx.hasSMSas(sms)
        #test getting efficiency with mass only
        self.assertAlmostEqual(tx.getEfficiencyFor(smsMatch), 0.12396,3)

        #test getting efficiency with SMS and reweighted efficiency
        c1.totalwidth = 1e-17*GeV
        self.assertAlmostEqual(tx.getEfficiencyFor(smsMatch), 0.08697,3)

    def testGetValueFor(self):

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM1.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        mass = [150.0,150.0]
        #test getting UL with mass only
        self.assertEqual(tx.txnameData.getValueFor(mass), 0.21496)

        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model( BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])
        sms = fromString("[[],[]]",finalState = ['C1+','C1-'],  model=model)
        c1m = model.getParticle(label='C1-')
        c1p = model.getParticle(label='C1+')
        c1p.mass = 150.0*GeV
        c1m.mass = 150.0*GeV
        c1p.totalwidth = 1e-17*GeV
        c1m.totalwidth = 1e-17*GeV
        smsMatch = tx.hasSMSas(sms)
        eff = tx.getEfficiencyFor(smsMatch)

        self.assertAlmostEqual(eff,0.49*0.21496,3)

    def testCoordinateTransf(self):
        """ test the transformation of data into coordinates, back into data """

        #Test with a regular mass array:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"],
                    datasetIDs=[None], txnames=["T2bb" ] )
        txname = expRes[0].datasets[0].txnameList[0] # T2bb
        initial = [300.,100.,300.,100.]
        coords = txname.txnameData.PCAtransf(initial)
        data = txname.txnameData.inversePCAtransf(coords)
        data = np.array(data)
        initial = np.array(initial)
        self.assertEqual(data.shape,initial.shape)
        diff = np.linalg.norm(data-initial)
        self.assertAlmostEqual(diff, 0.,2)

        #Test with with a mass array containing tuples:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2016-08"],
                    datasetIDs=[None], txnames=["T5Disp" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [300.,100.,300.,100.,1e-16,1e-16]
        coords = txname.txnameData.PCAtransf(initial)
        data = txname.txnameData.inversePCAtransf(coords)
        data = np.array(data)
        initial = np.array(initial)
        self.assertEqual(data.shape,initial.shape)
        diff = np.linalg.norm(data-initial)
        self.assertAlmostEqual(diff, 0.,2)


    def testCoordinateTransfInclusive(self):
        """ test the transformation of data into coordinates, back into data """

        #Test with an inclusive result:
        expRes = database.getExpResults(analysisIDs=["CMS-EXO-13-006"],
                    datasetIDs=['c000'], txnames=["THSCPM4" ])
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [300.,100.]
        coords = txname.txnameData.PCAtransf(initial)
        data = txname.txnameData.inversePCAtransf(coords)
        data = np.array(data)
        initial = np.array(initial)
        self.assertEqual(data.shape,initial.shape)
        diff = np.linalg.norm(data-initial)
        self.assertAlmostEqual(diff, 0.,2)

        #Test with another type of inclusive result:
        expRes = database.getExpResults(analysisIDs=["CMS-EXO-13-006"],
                    datasetIDs=['c000'], txnames=["THSCPM7" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [300.,50.,200.,100.,1.]
        coords = txname.txnameData.PCAtransf(initial)
        data = txname.txnameData.inversePCAtransf(coords)
        data = np.array(data)
        initial = np.array(initial)
        self.assertEqual(data.shape,initial.shape)
        diff = np.linalg.norm(data-initial)
        self.assertAlmostEqual(diff, 0.,2)


if __name__ == "__main__":
    unittest.main()
