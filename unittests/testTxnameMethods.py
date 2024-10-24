#!/usr/bin/env python3

"""
.. module:: testTxNameMethods
   :synopsis: Test some of the methods of the TxName class

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
import numpy as np
import unum
from unitTestHelpers import getNodesIndices, getEdges

filePath = './database/13TeV/CMS/CMS-PAS-SUS-15-002/data/T1.txt'
globalInfo = Info('./database/13TeV/CMS/CMS-PAS-SUS-15-002/globalInfo.txt')
infoObj = Info('./database/13TeV/CMS/CMS-PAS-SUS-15-002/data/dataInfo.txt')
databaseParticles = finalStates
tx = TxName(filePath,globalInfo,infoObj,databaseParticles)

class TestTxnameMethods(unittest.TestCase):
    
    
    def test_txload(self):

        sms = list(tx.smsMap.keys())[0]
        nodes_and_indices = getNodesIndices(sms)
        edges = getEdges(sms)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('MET', 3), 
                                          ('jet', 4), ('jet', 5), 
                                          ('MET', 6), ('jet', 7), 
                                          ('jet', 8)])
        self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'MET'), 
                              ('anyBSM', 'jet'), ('anyBSM', 'jet'), 
                              ('anyBSM', 'jet'), ('anyBSM', 'jet')])
        
        dataMap = {0 : (1, 'mass', 1.00E+00*GeV),
                   1 : (3, 'mass', 1.00E+00*GeV),
                   2 : (2, 'mass', 1.00E+00*GeV),
                   3 : (6, 'mass', 1.00E+00*GeV)}        
        self.assertEqual(tx.dataMap,dataMap)

        arrayMap = {0 : ((0, 0, 0), 'mass', 1.00E+00*GeV, 1),
                    1 : ((0, 1, 0), 'mass', 1.00E+00*GeV, 3),
                    2 : ((1, 0, 0), 'mass', 1.00E+00*GeV, 2),
                    3 : ((1, 1, 0), 'mass', 1.00E+00*GeV, 6)}
        self.assertEqual(tx._arrayMap,arrayMap)
        self.assertEqual(tx.y_unit,fb)

    def test_transform_point(self):

        massPoint = [[100*GeV,50*GeV],[0.2*TeV,10*GeV]]
        self.assertEqual(tx.transformPoint(massPoint),
                         [100.0, 50.0, 200.0, 10.0])
        
        data = [ [[[100*GeV,50*GeV],[200*GeV,10*GeV]],0.1*fb],  [[[400*GeV,150*GeV],[300*GeV,30*GeV]],10.*pb]]
        xvalues,yvalues = tx.transformData(data)
        xyvalues = list(zip([list(x) for x in xvalues],yvalues))
        xy = [([100.,50.,200.,10.],0.1),([400.,150.,300.,30.], 10000.0)]
        self.assertEqual(xyvalues,xy)

    def test_from_model(self):

        slhafile="./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)
        expSMS = ExpSMS.from_string("[[['u','u~']],[['d','d~']]]",model=model,
                    intermediateState=[['gluino'],['gluino']],finalState=['N1','N2'])
        # Hack to create a theory element from a string:
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        sms = list(tx.smsMap.keys())[0]
        smsMatch = sms.matchesTo(treeA)

        nodes_and_indices = getNodesIndices(smsMatch)
        edges = getEdges(smsMatch)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                          ('gluino', 2), ('N1', 3), 
                                          ('q', 4), ('q', 5), 
                                          ('N2', 6), ('q', 7), 
                                          ('q', 8)])
        self.assertEqual(edges,[('PV', 'gluino'), ('PV', 'gluino'), 
                              ('gluino', 'N1'), ('gluino', 'N2'), 
                              ('gluino', 'q'), ('gluino', 'q'), 
                              ('gluino', 'q'), ('gluino', 'q')])
        
        masses = [None, 5.77E+02*GeV, 5.77E+02*GeV, 
                  6.81E+01*GeV, 0.00E+00*MeV, 0.00E+00*MeV, 
                  1.353E+02*GeV, 0.00E+00*MeV, 0.00E+00*MeV]
        for im,m in enumerate(smsMatch.mass):
            if m is None:
                self.assertIs(masses[im],None)
                continue
            self.assertAlmostEqual(m.asNumber(GeV),
                             masses[im].asNumber(GeV),places=1)
            

        smsData = tx.getDataFromSMS(smsMatch)
        self.assertEqual(smsData,[577.0, 68.1, 577.0, 135.3])
        v = tx.txnameData.getValueFor(smsData)
        self.assertIs(v,None)
        reweigtF = tx.getReweightingFor(smsMatch)
        self.assertIs(reweigtF,None)

    def test_reweight(self):

        slhafile="./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)
        expSMS = ExpSMS.from_string("[[['u','u~']],[['d','d~']]]",
                                    model=model,             intermediateState=[['gluino'],['gluino']],finalState=['N1','N1'])
        # Hack to create a theory element from a string:
        smsB = TheorySMS()
        smsB.add_nodes_from(expSMS.nodes)
        smsB.add_edges_from(expSMS.edgeIndices)
        smsB.prodXSec = 1.0*fb
        smsB.maxWeight = 1.0*fb
        smsB.setGlobalProperties()

        sms = list(tx.smsMap.keys())[0]
        smsMatch = sms.matchesTo(smsB)

        nodes_and_indices = getNodesIndices(smsMatch)
        edges = getEdges(smsMatch)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                          ('gluino', 2), ('N1', 3), 
                                          ('q', 4), ('q', 5), 
                                          ('N1', 6), ('q', 7), 
                                          ('q', 8)])
        self.assertEqual(edges,[('PV', 'gluino'), ('PV', 'gluino'), 
                              ('gluino', 'N1'), ('gluino', 'N1'), 
                              ('gluino', 'q'), ('gluino', 'q'), 
                              ('gluino', 'q'), ('gluino', 'q')])

        masses = [None, 5.77E+02*GeV, 5.77E+02*GeV, 
                  6.81E+01*GeV, 0.00E+00*MeV, 
                  0.00E+00*MeV, 6.81E+01*GeV, 
                  0.00E+00*MeV, 0.00E+00*MeV]

        for im,m in enumerate(smsMatch.mass):
            if m is None:
                self.assertIs(masses[im],None)
                continue
            self.assertAlmostEqual(m.asNumber(GeV),
                             masses[im].asNumber(GeV),places=1)
            
        smsData = tx.getDataFromSMS(smsMatch)
        self.assertEqual(smsData,[577.0, 68.1, 577.0, 68.1])
        v = tx.txnameData.getValueFor(smsData)
        self.assertIs(v,None)
        reweigtF = tx.getReweightingFor(smsMatch)
        self.assertEqual(reweigtF,1.0)

        smsData = [650.,203.,650.,203.]
        v = tx.txnameData.getValueFor(smsData)
        self.assertAlmostEqual(v,195.8,places=1)
        # Change gluino and neutralino masses to get a result:
        gluino = model.getParticlesWith(label='gluino')[0]
        N1 = model.getParticlesWith(label='N1')[0]
        gluino.mass = 650*GeV
        N1.mass = 203*GeV
        v = tx.getULFor(smsMatch)
        self.assertAlmostEqual(v.asNumber(fb),195.8,places=1)

        p = tx.txnameData.PCAtransf(smsData)
        # print(p)
        pInv = tx.txnameData.inversePCAtransf(p)
        pInv = [tx.txnameData.round_to_n(x,5) for x in pInv[:]]
        self.assertEqual(pInv,[650.0, 203.0, 650.0, 203.0])
        self.assertEqual(pInv,smsData)

        pInv = tx.inverseTransformPoint(smsData)
        masses = [[6.50E+02*GeV, 2.03E+02*GeV], [6.50E+02*GeV, 2.03E+02*GeV]]
        for ibr,br in enumerate(pInv):
            for im,m in enumerate(br):
                self.assertAlmostEqual(m.asNumber(GeV),
                                masses[ibr][im].asNumber(GeV),places=1)
                
        # Check reweighting
        # Change gluino and neutralino widths:
        gluino.totalwidth = 1e-13*GeV
        N1.totalwidth = 1e-20*GeV
        Leff_inner = 0.000769
        Leff_outer = 7.0
        hc = 197.327*1e-18
        Flong = (np.exp(-N1.totalwidth.asNumber(GeV)*Leff_outer/hc))**2
        Fprompt = (1. - np.exp(-gluino.totalwidth.asNumber(GeV)*Leff_inner/hc))**2
        rFexp = 1/(Flong*Fprompt)
        rF = tx.getReweightingFor(smsMatch)
        self.assertAlmostEqual(rF,rFexp,places=4)
        ul = tx.getULFor(smsMatch)
        self.assertAlmostEqual(ul.asNumber(fb),196*9.6,places=0)


if __name__ == "__main__":
    unittest.main()                         