#!/usr/bin/env python3

"""
.. module:: testSMSComb
   :synopsis: Testing comparison of SMS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.experiment.defaultFinalStates import finalStates
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV
from unitTestHelpers import getNodesIndices, getEdges

slhafile="./testFiles/slha/lightEWinos.slha"
model = Model(BSMList,SMList)
model.updateParticles(inputFile=slhafile, ignorePromptQNumbers=['spin'])


class TestSMSComparison(unittest.TestCase):

    
    
    def test_sms_comp(self):


      smsA = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > W+,higgs), (H-(4) > W-,higgs), (H-(2) > W-,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(smsA, model=model)
      treeA = TheorySMS()
      treeA.add_nodes_from(expSMS.nodes)
      treeA.add_edges_from(expSMS.edgeIndices)
      treeA.prodXSec = 1.0*fb
      treeA.maxWeight = 1.0*fb
      treeA.setGlobalProperties()

      smsB = "(PV > H0(1),H-(2)), (H-(2) > W-,higgs), (H0(1) > H-(3),H+(4)), (H-(3) > higgs,W-), (H+(4) > W+,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(smsB, model=model)
      treeB = TheorySMS()
      treeB.add_nodes_from(expSMS.nodes)
      treeB.add_edges_from(expSMS.edgeIndices)
      treeB.prodXSec = 1.0*fb
      treeB.maxWeight = 1.0*fb
      treeB.setGlobalProperties()

      self.assertTrue(treeA == treeB)

      smsB = "(PV > H0(1),H-(2)), (H-(2) > W-,higgs), (H0(1) > H-(3),H+(4)), (H-(3) > higgs,W-), (H+(4) > W+,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(smsB, model=model)
      matched = expSMS.matchesTo(treeA)

      nodes_and_indices = getNodesIndices(matched)
      edges = getEdges(matched)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('H0', 1), 
                                          ('H-', 2), ('H-', 3), ('H+', 4),
                                          ('W-', 5), ('higgs', 6),
                                          ('higgs', 7), ('W-', 8),
                                          ('W+', 9), ('higgs', 10)])
      self.assertEqual(edges,[('H+', 'W+'), ('H+', 'higgs'),
                              ('H-', 'W-'), ('H-', 'W-'), 
                              ('H-', 'higgs'), ('H-', 'higgs'),
                              ('H0', 'H+'), ('H0', 'H-'), 
                              ('PV', 'H-'), ('PV', 'H0')])

            
      smsC = "(PV > H0(1),H-(2)), (H-(2) > W+,higgs), (H0(1) > H-(3),H+(4)), (H-(3) > higgs,W-), (H+(4) > W+,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(smsC, model=model)
      treeC = TheorySMS()
      treeC.add_nodes_from(expSMS.nodes)
      treeC.add_edges_from(expSMS.edgeIndices)
      treeC.prodXSec = 1.0*fb
      treeC.maxWeight = 1.0*fb
      treeC.setGlobalProperties()

      self.assertFalse(treeC == treeA)


    def test_comp2(self):

      sms = "(PV > gluino(1),gluino(2)), (gluino(1) > N1,b,b), (gluino(2) > g,N3(3)), (N3(3) > N1,Z)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(sms, model=model)
      treeA = TheorySMS()
      treeA.add_nodes_from(expSMS.nodes)
      treeA.add_edges_from(expSMS.edgeIndices)
      treeA.prodXSec = 1.0*fb
      treeA.maxWeight = 1.0*fb
      treeA.setGlobalProperties()

      sms = "(PV > gluino(1),gluino(2)), (gluino(1) > g,N3(3)), (gluino(2) > N1,b,b), (N3(3) > N1,Z)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(sms, model=model)
      treeB = TheorySMS()
      treeB.add_nodes_from(expSMS.nodes)
      treeB.add_edges_from(expSMS.edgeIndices)
      treeB.prodXSec = 1.0*fb
      treeB.maxWeight = 1.0*fb
      treeB.setGlobalProperties()

      self.assertTrue(treeA == treeB)

    def test_fromStr(self):
       
      sms1 = ExpSMS.from_string("(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),jet), (anyBSM(2) > anyBSM(4),jet), (anyBSM(3) > MET,Z), (anyBSM(4) > MET,W)",                            
             model=finalStates)
      sms2 = ExpSMS.from_string("(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),jet), (anyBSM(2) > anyBSM(4),jet), (anyBSM(3) > MET,W), (anyBSM(4) > MET,Z)",
             model=finalStates)
       
      self.assertTrue(sms1 == sms2)

      sms1 = ExpSMS.from_string("(PV > C1(1),MET), (C1(1) > anySM,MET)",model=finalStates)
      sms2 = ExpSMS.from_string("(PV > Hpm(1),MET), (Hpm(1) > anySM,MET)",model=finalStates)

      self.assertFalse(sms1 == sms2)

if __name__ == "__main__":
    unittest.main()                         
