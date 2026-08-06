#!/usr/bin/env python3

"""
.. module:: testTxNameConstraint
   :synopsis: Test the loading of txname constraints

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

filePath = './database/8TeV/ATLAS/ATLAS-SUSY-2013-12/data/TChiChipmSlepL.txt'
globalInfo = Info('./database/8TeV/ATLAS/ATLAS-SUSY-2013-12/globalInfo.txt')
infoObj = Info('./database/8TeV/ATLAS/ATLAS-SUSY-2013-12/data/dataInfo.txt')
databaseParticles = finalStates
tx = TxName(filePath,globalInfo,infoObj,databaseParticles)

class TestTxnameConstraint(unittest.TestCase):
    
    
    def test_old_constraint(self):

        self.assertEqual(tx._constraint,"2.*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])")
        self.assertEqual(tx._condition,["Csim([[['L'],['L']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]])", " Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['ta'],['ta']],[['nu'],['L']]])", " Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['ta'],['ta']],[['L'],['nu']]])", "Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['ta']]])", " Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['ta'],['nu']]])", "Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['e'],['e']],[['nu'],['L']]])", " Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['e'],['e']],[['L'],['nu']]])", " Cgtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['e']]])", " Cgtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['e'],['nu']]])"])

        self.assertEqual(tx._constraintFunc,"2.*(sms_1+sms_2)")
        self.assertEqual(len(tx._conditionsList),9)

        smsDicts = [{'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_1', 
                      '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,L)': 'sms_2'},                       
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,L)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),ta), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,ta), (anyBSM(7) > MET,L)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),ta), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,ta), (anyBSM(7) > MET,nu)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,L)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,ta)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),ta), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,L)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),e), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,e), (anyBSM(7) > MET,L)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),e), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,e), (anyBSM(7) > MET,nu)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,L)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),nu), (anyBSM(3) > MET,L), (anyBSM(7) > MET,e)': 'sms_2'}, 
                    {'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),L), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_1', 
                     '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),L), (anyBSM(2) > anyBSM(7),e), (anyBSM(3) > MET,L), (anyBSM(7) > MET,nu)': 'sms_2'}]
        condList = ['Csim(sms_1,sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)','Cgtr(sms_1,3.*sms_2)']
        for ic,cond in enumerate(tx._conditionsList):
            self.assertEqual(len(cond),1)
            c = list(cond.keys())[0]
            cDict = list(cond.values())[0]
            cDict = {str(k) : v for k,v in list(cDict.items())}
            self.assertEqual(c,condList[ic])
            self.assertEqual(cDict,smsDicts[ic])


    def test_smsMap(self):

        smsDict = {'sms_1' : 
                   {'nodes' : [('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('anyBSM', 3), 
                                          ('L', 4), ('MET', 5), 
                                          ('L', 6), ('anyBSM', 7), 
                                          ('L', 8), ('MET', 9), 
                                          ('nu', 10)],
                    'edges' : [('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'L'), ('anyBSM', 'L'), 
                              ('anyBSM', 'L'), ('anyBSM', 'MET'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'anyBSM'), 
                              ('anyBSM', 'anyBSM'), ('anyBSM', 'nu')]},
                   'sms_2':                   
                   {'nodes' : [('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('anyBSM', 3), 
                                          ('L', 4), ('MET', 5), 
                                          ('L', 6), ('anyBSM', 7), 
                                          ('nu', 8), ('MET', 9), 
                                          ('L', 10)],
                    'edges' : [('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'L'), ('anyBSM', 'L'), 
                              ('anyBSM', 'L'), ('anyBSM', 'MET'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'anyBSM'), 
                              ('anyBSM', 'anyBSM'), ('anyBSM', 'nu')]}}

        self.assertEqual(len(tx.smsMap),2)
        self.assertEqual(sorted(tx.smsMap.values()),['sms_1','sms_2'])
        for sms,smsLabel in tx.smsMap.items():
            nodes_and_indices = getNodesIndices(sms)
            edges = getEdges(sms)
            self.assertEqual(nodes_and_indices,smsDict[smsLabel]['nodes'])
            self.assertEqual(edges,smsDict[smsLabel]['edges'])


    def test_macthing(self):

        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)

        expSMS = ExpSMS.from_string('(PV > C1+(1),N2(2)), (C1+(1) > nu,se_L(3)), (N2(2) > e+,sne_L(4)), (se_L(3) > e-,N1), (sne_L(4) > e-,N1)',
             model=model)
        # Hack to create a theory element from a string:
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        expSMS = ExpSMS.from_string('(PV > C1+(1),N2(2)), (C1+(1) > e+,sne_L(3)), (N2(2) > e+,sne_L(4)), (sne_L(3) > nue,N1), (sne_L(4) > e-,N1)',
             model=model)
        # Hack to create a theory element from a string:
        treeB = TheorySMS()
        treeB.add_nodes_from(expSMS.nodes)
        treeB.add_edges_from(expSMS.edgeIndices)
        treeB.prodXSec = 1.0*fb
        treeB.maxWeight = 1.0*fb
        treeB.setGlobalProperties()


        smsMatchA = tx.hasSMSas(treeA)
        smsMatchB = tx.hasSMSas(treeB)
        smsMatchA.weight = 100*fb
        smsMatchB.weight = 5*pb
        


        self.assertEqual(smsMatchA.txlabel,'sms_2')
        nodes_and_indices = getNodesIndices(smsMatchA)
        edges = getEdges(smsMatchA)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('N2', 1), 
                                          ('C1+', 2), ('sne_L', 3), 
                                          ('e+', 4), ('N1', 5), 
                                          ('e-', 6), ('se_L', 7), 
                                          ('nu', 8), ('N1', 9), 
                                          ('e-', 10)])
        self.assertEqual(edges,[('C1+', 'nu'), ('C1+', 'se_L'), 
                              ('N2', 'e+'), ('N2', 'sne_L'), 
                              ('PV', 'C1+'), ('PV', 'N2'), 
                              ('se_L', 'N1'), ('se_L', 'e-'), 
                              ('sne_L', 'N1'), ('sne_L', 'e-')])
        
        
        self.assertEqual(smsMatchB.txlabel,'sms_1')
        nodes_and_indices = getNodesIndices(smsMatchB)
        edges = getEdges(smsMatchB)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('N2', 1), 
                                          ('C1+', 2), ('sne_L', 3), 
                                          ('e+', 4), ('N1', 5), 
                                          ('e-', 6), ('sne_L', 7), 
                                          ('e+', 8), ('N1', 9), 
                                          ('nu', 10)])
        self.assertEqual(edges,[('C1+', 'e+'), ('C1+', 'sne_L'), 
                              ('N2', 'e+'), ('N2', 'sne_L'), 
                              ('PV', 'C1+'), ('PV', 'N2'), 
                              ('sne_L', 'N1'), ('sne_L', 'N1'), 
                              ('sne_L', 'e-'), ('sne_L', 'nu')])

        allSMS = [smsMatchA,smsMatchB]
        constraint = tx.evalConstraintFor(allSMS)
        self.assertAlmostEqual(constraint.asNumber(pb),10.2,places=1)
        self.assertTrue(constraint == 2*(smsMatchA.weight+smsMatchB.weight))
   
        conds = [0.9607, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5]
        for ic,c in enumerate(tx.evalConditionsFor(allSMS)):
            self.assertAlmostEqual(c,conds[ic],places=3)

    def test_axes_conversion(self):


        # Test a trivial example
        axesMap = [{0: 'x', 1: '0.5*x+0.5*y', 2: 'y', 3: 'x', 4 : '0.5*x+0.5*y', 5 : 'y'}]
        self.assertTrue(isinstance(tx.axesMap,list))
        self.assertEqual(len(tx.axesMap),1)
        for iax, ax in enumerate(sorted(tx.axesMap, key = lambda x: str(x))):
            self.assertEqual(ax,axesMap[iax])

        # Test an example with widths
        filePath = './database_extra/13TeV/CMS/CMS-EXO-19-010-eff/SR_nlay5/TDTM1F.txt'
        globalInfo = Info('./database_extra/13TeV/CMS/CMS-EXO-19-010-eff/globalInfo.txt')
        infoObj = Info('./database_extra/13TeV/CMS/CMS-EXO-19-010-eff/SR_nlay5/dataInfo.txt')
        databaseParticles = finalStates
        txB = TxName(filePath,globalInfo,infoObj,databaseParticles)

        axesMap = sorted([{0: 'x', 1: 'x-1.5', 2: 'x', 3: 'x-1.5', 4: 'y', 5: 'y'}, 
                   {0: 'x', 1: 'x', 2: 'x', 3: 'x', 4: 'y', 5: 'y'}], key = lambda x: str(x))
        self.assertTrue(isinstance(txB.axesMap,list))
        self.assertEqual(len(txB.axesMap),2)
        for iax, ax in enumerate(sorted(txB.axesMap, key = lambda x: str(x))):
            self.assertEqual(ax,axesMap[iax])

        # Test a non-trivial case
        filePath = './database_extra/13TeV/ATLAS/ATLAS-SUSY-2016-32-eff-trim/SR1FULL_175/THSCPM4.txt'
        globalInfo = Info('./database_extra/13TeV/ATLAS/ATLAS-SUSY-2016-32-eff-trim/globalInfo.txt')
        infoObj = Info('./database_extra/13TeV/ATLAS/ATLAS-SUSY-2016-32-eff-trim/SR1FULL_175/dataInfo.txt')
        databaseParticles = finalStates
        txC = TxName(filePath,globalInfo,infoObj,databaseParticles)

        axesMap = [{0: 'x', 1: 'y', 2: 'w'}]
        self.assertTrue(isinstance(txC.axesMap,list))
        self.assertEqual(len(txC.axesMap),1)
        for iax, ax in enumerate(sorted(txC.axesMap, key = lambda x: str(x))):
            self.assertEqual(ax,axesMap[iax])


if __name__ == "__main__":
    unittest.main()                         