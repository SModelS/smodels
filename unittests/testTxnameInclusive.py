#!/usr/bin/env python3

"""
.. module:: testTxNameInclusive
   :synopsis: Test the use of inclusive topologies

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
from smodels.base.physicsUnits import fb, GeV, TeV, pb
from smodels.experiment.txnameObj import TxName, TxNameData
from smodels.experiment.infoObj import Info
from unitTestHelpers import getNodesIndices, getEdges

filePath = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM4.txt'
globalInfo = Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
infoObj = Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/dataInfo.txt')
databaseParticles = finalStates
tx = TxName(filePath,globalInfo,infoObj,databaseParticles)

class TestTxnameInclusive(unittest.TestCase):
    
    
    def test_load(self):
        sms = list(tx.smsMap.keys())[0]
        nodes_and_indices = getNodesIndices(sms)
        edges = getEdges(sms)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('Inclusive', 1), 
                                          ('anyBSM', 2), ('MET', 3), 
                                          ('*anySM', 4), ('HSCP', 5), 
                                          ('anySM', 6)])
        self.assertEqual(edges,[('Inclusive', '*anySM'), 
                                ('Inclusive', 'MET'), 
                              ('PV', 'Inclusive'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'HSCP'), ('anyBSM', 'anySM')])
        
        dataMap = {0 : (2,'mass',GeV), 1 : (5, 'mass', GeV)}
        self.assertEqual(tx.dataMap,dataMap)

        arrayMap = {0 : ((1, 0, 0), 'mass', GeV, 2), 
                    1 : ((1, 1, 0), 'mass', GeV, 5)}
        self.assertEqual(tx._arrayMap,arrayMap)

        self.assertEqual(tx.y_unit,1.0)

    def test_transform(self):

        massPoint = [[(100*GeV,1e-3*GeV),50*GeV],[(0.2*TeV,1e6*GeV),10*GeV]]
        self.assertEqual(tx.transformPoint(massPoint),[200.0, 10.0])
        massPoint = ['*',[(0.1*TeV,1e6*GeV),20*GeV]]
        self.assertEqual(tx.transformPoint(massPoint),[100.0, 20.0])

        data = [ [['*',[200*GeV,10*GeV]],0.1*fb],  
        [['*',[300*GeV,30*GeV]],10.*pb]]
        xvalues,yvalues = tx.transformData(data)
        xyvalues = list(zip([list(x) for x in xvalues],yvalues))
        xy = [([200., 10.],0.1),([300.,30.], 10000.0)]
        self.assertEqual(xyvalues,xy)   

    def test_matching(self):

        sms = list(tx.smsMap.keys())[0]

        slhafile="./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)
        # Hack to create a theory smsement from a string:
        expSMS = ExpSMS.from_string("[[['W-']],[['d','d~']]]",model=model,
                    intermediateState=[['N3'],['N2']],finalState=['C1+','N1'])
        smsB = TheorySMS()
        smsB.add_nodes_from(expSMS.nodes)
        smsB.add_edges_from(expSMS.edgeIndices)
        smsB.prodXSec = 1.0*fb
        smsB.maxWeight = 1.0*fb
        smsB.setGlobalProperties()
        smsMatch = sms.matchesTo(smsB)

        nodes_and_indices = getNodesIndices(smsMatch)
        edges = getEdges(smsMatch)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('N2', 1), 
                                          ('N3', 2), ('N1', 3), 
                                          ('q', 4), ('C1+', 5), 
                                          ('W-', 6), ('q', 7)])
        self.assertEqual(edges,[('N2', 'N1'), ('N2', 'q'), 
                              ('N2', 'q'), ('N3', 'C1+'), 
                              ('N3', 'W-'), ('PV', 'N2'), 
                              ('PV', 'N3')])
        
        smsData = tx.getDataFromSMS(smsMatch)
        self.assertEqual(smsData,[266.5, 134.4])
        v = tx.txnameData.getValueFor(smsData)
        self.assertAlmostEqual(1.43e-1,v,places=3)
        reweigtF = tx.getReweightingFor(smsMatch)
        self.assertEqual(reweigtF,0.0)

        # Change chargino width to get a result:
        C1p = model.getParticle(label='C1+')
        C1m = model.getParticle(label='C1-')
        C1p.totalwidth = 5e-17*GeV
        C1m.totalwidth = 5e-17*GeV
        v = tx.getEfficiencyFor(smsMatch)
        self.assertAlmostEqual(2.427e-02,v,places=5)

        p = tx.txnameData.PCAtransf(smsData)
        pInv = tx.txnameData.inversePCAtransf(p)
        pInv = [tx.txnameData.round_to_n(x,5) for x in pInv[:]]
        pInvB = [tx.txnameData.round_to_n(x,5) for x in smsData]
        self.assertEqual(smsData,[266.5, 134.4])
        self.assertEqual(pInv, pInvB)

        m = tx.inverseTransformPoint(smsData)
        self.assertEqual(len(m),2)
        self.assertEqual(m[0],[])
        massV = [266.5,134.4]
        for imass,mass in enumerate(m[1]):
            self.assertAlmostEqual(mass.asNumber(GeV),massV[imass],places=1)

if __name__ == "__main__":
    unittest.main()                         