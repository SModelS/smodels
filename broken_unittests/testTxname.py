#!/usr/bin/env python3

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from smodels.share.models import mssm
from smodels.theory.element import Element
from smodels.experiment import infoObj
from smodels.experiment.txnameObj import TxName, rescaleWidth, unscaleWidth
from smodels.tools.physicsUnits import GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import flattenArray
from databaseLoader import database
from smodels.experiment.defaultFinalStates import finalStates
import numpy as np
import unittest

c1 = mssm.c1
n1 = mssm.n1

class TxTest(unittest.TestCase):
    def testWidthTransformation(self):
        a = 1e-10*GeV
        b = unscaleWidth(rescaleWidth(a))
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 1e-50*GeV
        b = unscaleWidth(rescaleWidth(a))
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = float('inf')*GeV
        b = unscaleWidth(rescaleWidth(a))
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 1e-25*GeV
        b = unscaleWidth(rescaleWidth(a))
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))
        a = 0.*GeV
        b = unscaleWidth(rescaleWidth(a))
        self.assertAlmostEqual(b.asNumber(GeV), a.asNumber(GeV))

    def testTxnameElements(self):

        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        el = Element(info="[[*],[]]",finalState = ['MET','HSCP'], model = finalStates)

        self.assertTrue(len(tx._topologyList.getElements()), 1)
        self.assertEqual(tx._topologyList.getElements()[0], el)

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

        errstr = "BSM particle ``non-MET'' has not been defined in databaseParticles.py"
        self.assertEqual(gotError.args[0], errstr)

    def testgetEffFor(self):

        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        self.assertFalse(hasattr(tx._topologyList.getElements()[0], 'mass'))

        el = Element(info="[[[e+]],[]]",finalState = ['HSCP','MET'], model = finalStates)
        setattr(n1, 'mass', 200*GeV)
        setattr(c1, 'mass', 150*GeV)
        el.branches[0].oddParticles = [n1]
        el.branches[1].oddParticles = [c1]

        #test getting efficiency with mass only
        self.assertEqual(tx.getEfficiencyFor(el.mass), 0.12396)

        #test getting efficiency with element and reweighted efficiency
        setattr(n1, 'totalwidth', 0.*GeV)
        setattr(c1, 'totalwidth', 10**(-17)*GeV)
        self.assertAlmostEqual(tx.getEfficiencyFor(el), 0.08697,3)

    def testGetValueFor(self):

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM1.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo,finalStates)

        el = Element(info="[[],[]]",finalState = ['HSCP','HSCP'], model = finalStates)
        setattr(c1, 'mass', 150*GeV)
        el.branches[0].oddParticles = [c1]
        el.branches[1].oddParticles = [c1]

        #test getting UL with mass only
        self.assertEqual(tx.txnameData.getValueFor(el.mass), 0.21496)

        #test getting UL with element and reweighted efficiency
        setattr(c1, 'totalwidth', 10**(-17)*GeV)
        self.assertAlmostEqual(tx.txnameData.getValueFor(el),0.49*0.21496,3)

    def testCoordinateTransf(self):
        """ test the transformation of data into coordinates, back into data """

        #Test with a regular mass array:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"],
                    datasetIDs=[None], txnames=["T2bb" ] )
        txname=expRes[0].datasets[0].txnameList[0] # T2bb
        initial = [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ]
        coords=txname.txnameData.dataToCoordinates(
                initial, txname.txnameData._V, txname.txnameData.delta_x )
        data = txname.txnameData.coordinatesToData( coords, txname.txnameData._V,
                  txname.txnameData.delta_x)
        data = np.array(data)
        initial = np.array(initial)
        self.assertEqual(data.shape,initial.shape)
        dataFlat = np.array([x.asNumber(GeV) for x in flattenArray(data)])
        initialFlat = np.array([x.asNumber(GeV) for x in flattenArray(initial)])
        diff = np.linalg.norm(dataFlat-initialFlat)
        self.assertAlmostEqual(diff, 0.)

        #Test with with a mass array containing tuples:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2016-08"],
                    datasetIDs=[None], txnames=["T5Disp" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [[(300.*GeV,1e-16*GeV),100.*GeV], [(300.*GeV,1e-16*GeV),100.*GeV]]
        coords=txname.txnameData.dataToCoordinates(
                initial, txname.txnameData._V, txname.txnameData.delta_x )
        data = txname.txnameData.coordinatesToData(coords, txname.txnameData._V,
                  txname.txnameData.delta_x)
        data = np.array(data,dtype=object)
        initial = np.array(initial,dtype=object)
        self.assertEqual(data.shape,initial.shape)
        dataFlat = np.array([x.asNumber(GeV) for x in flattenArray(data)])
        initialFlat = np.array([x.asNumber(GeV) for x in flattenArray(initial)])
        diff = np.linalg.norm(dataFlat-initialFlat)
        self.assertAlmostEqual(diff, 0.)

    def testCoordinateTransfInclusive(self):
        """ test the transformation of data into coordinates, back into data """

        #Test with an inclusive result:
        expRes = database.getExpResults(analysisIDs=["CMS-EXO-13-006"],
                    datasetIDs=['c000'], txnames=["THSCPM4" ])
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [[(300.*GeV,1e-16*GeV),100.*GeV], [(300.*GeV,1e-16*GeV),100.*GeV]]
        coords=txname.txnameData.dataToCoordinates(
                initial, txname.txnameData._V, txname.txnameData.delta_x)
        data = txname.txnameData.coordinatesToData(coords, txname.txnameData._V,
                  txname.txnameData.delta_x)
        data = np.array(data,dtype=object)
        newInitial = ['*',[300.*GeV,100.*GeV]]
        newInitial = np.array(newInitial,dtype=object)
        self.assertEqual(data.shape,newInitial.shape)
        dataFlat = np.array([x.asNumber(GeV) if str(x) != '*' else -1 for x in flattenArray(data)])
        initialFlat = np.array([x.asNumber(GeV) if str(x) != '*' else -1 for x in flattenArray(newInitial)])
        diff = np.linalg.norm(dataFlat-initialFlat)
        self.assertAlmostEqual(diff, 0.)

        #Test with another type of inclusive result:
        expRes = database.getExpResults(analysisIDs=["CMS-EXO-13-006"],
                    datasetIDs=['c000'], txnames=["THSCPM7" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        initial = [[(300.*GeV,1e-16*GeV),50.*GeV], [(200.*GeV,1e-18*GeV),100.*GeV,1.*GeV]]
        coords=txname.txnameData.dataToCoordinates(
                initial, txname.txnameData._V, txname.txnameData.delta_x)
        data = txname.txnameData.coordinatesToData(coords, txname.txnameData._V,
                  txname.txnameData.delta_x)
        data = np.array(data,dtype=object)
        newInitial = [[300.*GeV,50.*GeV], [200.*GeV,100.*GeV,1.*GeV]]
        newInitial = np.array(newInitial,dtype=object)
        self.assertEqual(data.shape,newInitial.shape)
        dataFlat = np.array([x.asNumber(GeV) for x in flattenArray(data)])
        initialFlat = np.array([x.asNumber(GeV) for x in flattenArray(newInitial)])
        diff = np.linalg.norm(dataFlat-initialFlat)
        self.assertAlmostEqual(diff, 0.)


if __name__ == "__main__":
    unittest.main()
