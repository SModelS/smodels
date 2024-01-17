#!/usr/bin/env python3

"""
.. module:: testTxNameWidths
   :synopsis: Test the use of widths in TxNames

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.defaultFinalStates import finalStates
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV, TeV, pb, MeV
from smodels.experiment.txnameObj import TxName
from smodels.experiment.infoObj import Info
from smodels.experiment.expAuxiliaryFuncs import rescaleWidth
import numpy as np
from unitTestHelpers import getNodesIndices, getEdges


filePath = './database/13TeV/ATLAS/ATLAS-SUSY-2016-08/data/T5Disp.txt'
globalInfo = Info('./database/13TeV/ATLAS/ATLAS-SUSY-2016-08/globalInfo.txt')
infoObj = Info('./database/13TeV/ATLAS/ATLAS-SUSY-2016-08/data/dataInfo.txt')
databaseParticles = finalStates
tx = TxName(filePath,globalInfo,infoObj,databaseParticles)
sms = list(tx.smsMap.keys())[0]

class TestTxnameWidths(unittest.TestCase):
    
    def test_txload(self):

        nodes_and_indices = getNodesIndices(sms)
        edges = getEdges(sms)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('MET', 3), 
                                          ('q', 4), ('q', 5), 
                                          ('MET', 6), ('q', 7), 
                                          ('q', 8)])
        self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'MET'), 
                              ('anyBSM', 'q'), ('anyBSM', 'q'), 
                              ('anyBSM', 'q'), ('anyBSM', 'q')])
        
        dataMap = {0: (1, 'mass', 1.00E+00*GeV), 
                   1: (3, 'mass', 1.00E+00*GeV), 
                   2: (2, 'mass', 1.00E+00*GeV), 
                   3: (6, 'mass', 1.00E+00*GeV), 
                   4: (1, 'totalwidth', 1.00E+00*GeV), 
                   5: (2, 'totalwidth', 1.00E+00*GeV)}
        self.assertEqual(tx.dataMap,dataMap)

        arrayMap = {0: ((0, 0, 0), 'mass', 1.00E+00*GeV, 1), 
                    1: ((0, 1, 0), 'mass', 1.00E+00*GeV, 3), 
                    2: ((1, 0, 0), 'mass', 1.00E+00*GeV, 2), 
                    3: ((1, 1, 0), 'mass', 1.00E+00*GeV, 6), 
                    4: ((0, 0, 1), 'totalwidth', 1.00E+00*GeV, 1), 
                    5: ((1, 0, 1), 'totalwidth', 1.00E+00*GeV, 2)}

        self.assertEqual(tx._arrayMap,arrayMap)
        self.assertEqual(tx.y_unit,fb)
        
    def test_transform(self):

        massPoint = [[(100*GeV,1e-3*GeV),50*GeV],[(0.2*TeV,1e6*GeV),10*GeV]]
        massTransf = [100.0, 50.0, 200.0, 10.0, 62.1698, 82.8931]
        for im,m in enumerate(tx.transformPoint(massPoint)):
            self.assertAlmostEqual(m,massTransf[im],places=2)

        data = [ [[[(100*GeV,1e-3*GeV),50*GeV],[(200*GeV,1e6*GeV),10*GeV]],0.1*fb],  
        [[[(400*GeV,1*GeV),150*GeV],[(300*GeV,10*GeV),30*GeV]],10.*pb]]
        xvalues,yvalues = tx.transformData(data)
        xyvalues = list(zip([list(x) for x in xvalues],yvalues))
        xy = [([100.0, 50.0, 200.0, 10.0, 62.1698, 82.8931],0.1),
              ([400.0, 150.0, 300.0, 30.0, 69.0776, 71.3801], 10000.0)]
        for ixy,xyV in enumerate(xyvalues):            
            self.assertAlmostEqual(xyV[1],xy[ixy][1],places=3)
            for im,m in enumerate(xyV[0]):
                self.assertAlmostEqual(m,xy[ixy][0][im],places=3)

    def test_match(self):

        slhafile="./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string("[[['u','u~']],[['d','d~']]]",model=model,
                    intermediateState=[['C1+'],['C1-']],finalState=['N1','N1'])
        smsB = TheorySMS()
        smsB.add_nodes_from(expSMS.nodes)
        smsB.add_edges_from(expSMS.edgeIndices)
        smsB.prodXSec = 1.0*fb
        smsB.maxWeight = 1.0*fb
        smsB.setGlobalProperties()

        smsMatch = sms.matchesTo(smsB)
        nodes_and_indices = getNodesIndices(smsMatch)
        edges = getEdges(smsMatch)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('C1+', 1), 
                                          ('C1-', 2), ('N1', 3), 
                                          ('q', 4), ('q', 5), 
                                          ('N1', 6), ('q', 7), 
                                          ('q', 8)])
        self.assertEqual(edges,[('C1+', 'N1'), ('C1+', 'q'), 
                              ('C1+', 'q'), ('C1-', 'N1'), 
                              ('C1-', 'q'), ('C1-', 'q'), 
                              ('PV', 'C1+'), ('PV', 'C1-')])
        masses = [None, 1.344E+02*GeV, 1.344E+02*GeV, 6.81E+01*GeV, 0.00E+00*MeV, 0.00E+00*MeV, 6.81E+01*GeV, 0.00E+00*MeV, 0.00E+00*MeV]
        for im,m in enumerate(smsMatch.mass):
            if m is None:
                self.assertIs(masses[im],None)
                continue
            self.assertAlmostEqual(m.asNumber(GeV),
                             masses[im].asNumber(GeV),places=1)


        smsData = [np.round(x,3) for x in tx.getDataFromSMS(smsMatch)]
        massPoint = [134.4, 68.1, 134.4, 68.1, 59.696, 59.696]
        for im,m in enumerate(smsData):
            self.assertAlmostEqual(m,massPoint[im],places=1)
        v = tx.txnameData.getValueFor(smsData)
        self.assertIs(v,None)
        reweightF = tx.getReweightingFor(smsMatch)
        self.assertEqual(reweightF,1.0)

        smsData = [1450.,100.,1450.,100.,rescaleWidth(5e-17),rescaleWidth(5e-17)]
        v = tx.txnameData.getValueFor(smsData)
        self.assertAlmostEqual(v,1.547,places=2)

        # Change chargino width to get a result:
        C1p = model.getParticlesWith(label='C1+')[0]
        C1m = model.getParticlesWith(label='C1-')[0]
        N1 = model.getParticlesWith(label='N1')[0]
        C1p.totalwidth = 5e-17*GeV
        C1m.totalwidth = 5e-17*GeV
        C1p.mass = 1450.*GeV
        C1m.mass = 1450.*GeV
        N1.mass = 100.*GeV
        v = tx.getULFor(smsMatch)
        self.assertAlmostEqual(v.asNumber(fb),1.547,places=2)

        p = tx.txnameData.PCAtransf(smsData)
        pInv = tx.txnameData.inversePCAtransf(p)
        pInv = [tx.txnameData.round_to_n(x,5) for x in pInv[:]]
        massPoint = [1450.0, 100.0, 1450.0, 100.0, 31.543, 31.543]
        for im,m in enumerate(pInv):
            self.assertAlmostEqual(m,massPoint[im],places=2)
        self.assertTrue(pInv == [tx.txnameData.round_to_n(x,5) for x in smsData])

        massInv = tx.inverseTransformPoint(smsData)
        massRes = [[(1.45E+03*GeV, 5.00E-17*GeV), 1.00E+02*GeV], [(1.45E+03*GeV, 5.00E-17*GeV), 1.00E+02*GeV]]
        for ibr,br in enumerate(massInv):
            for im,m in enumerate(br):
                self.assertEqual(type(m),type(massRes[ibr][im]))
                if isinstance(m,tuple):
                    m1,w1 = m
                    m2,w2 = massRes[ibr][im]
                    self.assertAlmostEqual(m1.asNumber(GeV),m2.asNumber(GeV),places=2)
                    self.assertAlmostEqual(w1.asNumber(GeV),w2.asNumber(GeV),places=2)
                else:
                    self.assertAlmostEqual(m.asNumber(GeV),
                                            massRes[ibr][im].asNumber(GeV),places=2)

if __name__ == "__main__":
    unittest.main()                         