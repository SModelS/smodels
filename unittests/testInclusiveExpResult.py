#!/usr/bin/env python3

"""
.. module:: testInclusiveExpResult
   :synopsis: Testing creation and comparison of Inclusive SMS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
import os
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.experiment.databaseObj import Database
from smodels.experiment.txnameObj import TxName,TxNameData
from smodels.experiment.expSMS import ExpSMS
from smodels.experiment.infoObj import Info
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV
from databaseLoader import database as db
from smodels.base.inclusiveObjects import InclusiveValue
from unitTestHelpers import getNodesIndices, getEdges
db_extra = Database('unittestextra')


class TestInclusiveExpRes(unittest.TestCase):

    
    
    def test_expsms_inclusive(self):

        exp = db.getExpResults(analysisIDs='CMS-PAS-EXO-16-036', txnames='THSCPM2')[0]
        tx = exp.getTxNames()[0]
        self.assertEqual(len(tx.smsMap),1)
        expsms = list(tx.smsMap.keys())[0]
        nodes_and_indices = getNodesIndices(expsms)
        self.assertEqual(nodes_and_indices,
                         [('PV', 0), ('Inclusive', 1), ('HSCP', 2), 
                           ('MET', 3), ('*anySM', 4)])
        self.assertTrue(isinstance(expsms.canonName,InclusiveValue))
  


        # Test matching
        smsA = '(PV > gluino(1),sta_1), (gluino(1) > u,u~,N2(3)), (N2(3) > W-,sta_1~(4)), (sta_1~(4) > W+,N1)'
        slhafile="./testFiles/slha/longLived.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)
        smsA = ExpSMS.from_string(smsA,model=model)

        matched = tx.hasSMSas(smsA)
        nodes_and_indices = getNodesIndices(matched)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                            ('sta_1', 2), ('N1', 3), 
                                            ('W+', 4), ('N2', 5), ('sta_1~', 6), ('q', 7), 
                                            ('q', 8), ('W-', 9)])
        

        smsB = '(PV > sta_1,gluino(1)), (gluino(1) > u,u~,N1)'
        smsB = ExpSMS.from_string(smsB,model=model)
        matched = tx.hasSMSas(smsB)
        nodes_and_indices = getNodesIndices(matched)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                            ('sta_1', 2), ('N1', 3), 
                                            ('q', 4), ('q', 5)])
        
        edges = getEdges(matched)
        self.assertEqual(edges,sorted([('PV', 'gluino'), ('PV', 'sta_1'), 
                                       ('gluino', 'q'), ('gluino', 'N1'),
                                        ('gluino', 'q')]))
        

        smsC = '(PV > sta_1,gluino(1)), (gluino(1) > u,u~,N2(2)), (N2(2) > Z, N1)'
        smsC = ExpSMS.from_string(smsC,model=model)
        matched = tx.hasSMSas(smsC)
        nodes_and_indices = getNodesIndices(matched)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                            ('sta_1', 2), ('N1', 3), 
                                            ('Z', 4), ('N2', 5), 
                                            ('q', 6), ('q', 7)])
        
        smsD = '(PV > sta_1,gluino(1)), (gluino(1) > u,u~,N2(2)), (N2(2) > W-, sta_1)'
        smsD = ExpSMS.from_string(smsD,model=model)
        self.assertIs(tx.hasSMSas(smsD),None)

        smsE = '(PV > sta_1,sta_1)'
        smsE = ExpSMS.from_string(smsE,model=model)
        self.assertIs(tx.hasSMSas(smsE),None)

    def test_inclusive_dt(self):

        exp = db_extra.getExpResults(analysisIDs='CMS-EXO-19-010', txnames='TDTM1F')[0]
        tx = exp.getTxNames()[0]
        self.assertEqual(len(tx.smsMap),1)
        expsms = list(tx.smsMap.keys())[0]
        self.assertTrue(isinstance(expsms.canonName,InclusiveValue))
        nodes_and_indices = getNodesIndices(expsms)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('C1', 1), 
                                            ('C1', 2), ('MET', 3), 
                                            ('*anySM', 4), ('MET', 5), 
                                            ('*anySM', 6)])
        
        smsA = '(PV > C1-(1),C1+(2)), (C1+(2) > u,N1), (C1-(1) > W-,N1)'
        slhafile="./testFiles/slha/longLived.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1*GeV)
        smsA = ExpSMS.from_string(smsA,model=model)
        matched = tx.hasSMSas(smsA)
        nodes_and_indices = getNodesIndices(matched)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('C1-', 1), 
                                            ('C1+', 2), ('N1', 3), 
                                            ('W-', 4), ('N1', 5), 
                                            ('q', 6)])
        edges = getEdges(matched)
        self.assertEqual(edges,sorted([('PV', 'C1+'), ('PV', 'C1-'), 
                                       ('C1+', 'q'), ('C1+', 'N1'), 
                                       ('C1-', 'W-'), ('C1-', 'N1')]))
        

        smsA = '(PV > C1-(1),C1+(2)), (C1+(2) > u,N1), (C1-(1) > q,N1,e-,nu)'
        smsA = ExpSMS.from_string(smsA,model=model)
        matched = tx.hasSMSas(smsA)
        nodes_and_indices = getNodesIndices(matched)
        self.assertEqual(nodes_and_indices,[('PV', 0), ('C1-', 1), 
                                            ('C1+', 2), ('N1', 3), 
                                            ('q', 4), ('N1', 5), 
                                            ('q', 6), ('nu', 7), 
                                            ('e-', 8)])
        edges = getEdges(matched)
        self.assertEqual(edges,sorted([('C1+', 'N1'), ('C1+', 'q'), 
                                       ('C1-', 'N1'), ('C1-', 'e-'), 
                                       ('C1-', 'nu'), ('C1-', 'q'), 
                                       ('PV', 'C1+'), ('PV', 'C1-')]))
        

        smsA = '(PV > C1+(1),C1-(2)), (C1+(1) > u,N1), (C1-(2) > N2,e-,nu,N1)'
        smsA = ExpSMS.from_string(smsA,model=model)
        matched = tx.hasSMSas(smsA)
        self.assertIs(matched,None)
        
            
if __name__ == "__main__":
    unittest.main()                         