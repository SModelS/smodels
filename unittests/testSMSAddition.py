#!/usr/bin/env python3

"""
.. module:: testSMSAdd
   :synopsis: Testing combination of SMS

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
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV
from unitTestHelpers import getNodesIndices, getEdges

slhafile="./testFiles/slha/lightEWinos.slha"
model = Model(BSMList,SMList)
model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)



class TestSMSAddition(unittest.TestCase):

    
    
    def test_sms_add(self):

      smsA = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > W+,higgs), (H-(4) > W-,higgs), (H-(2) > W-,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(smsA, model=model)
      treeA = TheorySMS()
      treeA.add_nodes_from(expSMS.nodes)
      treeA.add_edges_from(expSMS.edgeIndices)
      treeA.prodXSec = 1.0*fb
      treeA.maxWeight = 1.0*fb
      treeA.setGlobalProperties()

      smsB = "(PV > gluino(1),gluino(2)), (gluino(2) > t+,t-), (gluino(1) > su_L(3),su_L~(4)), (su_L(3) > higgs,W+), (su_L~(4) > higgs,W-)"
      expSMS = ExpSMS.from_string(smsB, model=model)
      treeB = TheorySMS()
      treeB.add_nodes_from(expSMS.nodes)
      treeB.add_edges_from(expSMS.edgeIndices)
      treeB.prodXSec = 5.0*fb
      treeB.maxWeight = 5.0*fb
      treeB.setGlobalProperties()

      newTree = treeA+treeB
      nodes_and_indices = getNodesIndices(newTree)
      edges = getEdges(newTree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('H-/gluino', 1), 
                                          ('H0/gluino', 2), ('higgs/t-', 3), 
                                          ('W-/t+', 4), ('H-/su_L~', 5), 
                                          ('H+/su_L', 6), ('higgs', 7), ('W-', 8), 
                                          ('higgs', 9), ('W+', 10)])
      self.assertEqual(edges,[('H+/su_L', 'W+'), ('H+/su_L', 'higgs'), 
                              ('H-/gluino', 'W-/t+'), 
                              ('H-/gluino', 'higgs/t-'), 
                              ('H-/su_L~', 'W-'), 
                              ('H-/su_L~', 'higgs'), 
                              ('H0/gluino', 'H+/su_L'), 
                              ('H0/gluino', 'H-/su_L~'), 
                              ('PV', 'H-/gluino'), ('PV', 'H0/gluino')])

            
if __name__ == "__main__":
    unittest.main()                         