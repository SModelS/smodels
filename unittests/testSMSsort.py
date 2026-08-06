#!/usr/bin/env python3

"""
.. module:: testSMSsort
   :synopsis: Testing sorting of SMS objects

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



class TestSMSSort(unittest.TestCase):
    
    
    def test_sorting(self):
      stringProc = "(PV > gluino(1),gluino(2)), (gluino(1) > N3(3),g), (N3(3) > Z,N1), (gluino(2) > b,N1,b)"
      slhafile="./testFiles/slha/higgsinoStop.slha"
      model = Model(BSMList,SMList)
      model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      tree = TheorySMS()
      tree.add_nodes_from(expSMS.nodes)
      tree.add_edges_from(expSMS.edgeIndices)

      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                          ('gluino', 2), ('N3', 3), 
                                          ('g', 4), ('Z', 5), 
                                          ('N1', 6), ('b', 7), 
                                          ('N1', 8), ('b', 9)])
      self.assertEqual(edges,[('N3', 'N1'), ('N3', 'Z'), 
                              ('PV', 'gluino'), ('PV', 'gluino'), 
                              ('gluino', 'N1'), ('gluino', 'N3'), 
                              ('gluino', 'b'), ('gluino', 'b'), 
                              ('gluino', 'g')])
      
      tree.sort()
      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                          ('gluino', 2), ('N1', 3), 
                                          ('b', 4), ('b', 5), 
                                          ('g', 6), ('N3', 7), 
                                          ('N1', 8), ('Z', 9)])
      self.assertEqual(edges,[('N3', 'N1'), ('N3', 'Z'), 
                              ('PV', 'gluino'), ('PV', 'gluino'), 
                              ('gluino', 'N1'), ('gluino', 'N3'), 
                              ('gluino', 'b'), ('gluino', 'b'), 
                              ('gluino', 'g')])
      
      nodesCanonNames = [('PV', 0, 11101010011011010000),
                        ('gluino', 1, 11010100),
                        ('gluino', 2, 1101101000),
                        ('N1', 3, 10),
                        ('b', 4, 10),
                        ('b', 5, 10),
                        ('g', 6, 10),
                        ('N3', 7, 110100),
                        ('N1', 8, 10),
                        ('Z', 9, 10)]
      
      nodeListA = []
      for node,nodeIndex in zip(tree.nodes,tree.nodeIndices):
          nodeListA.append((str(node),nodeIndex,
                          tree.nodeCanonName(nodeIndex)))
      self.assertEqual(nodesCanonNames,nodeListA)


      stringProc = "(PV > gluino(1),gluino(2)), (gluino(1) > N1,b,b), (gluino(2) > N3(3),g), (N3(3) > N1,Z)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      tree = TheorySMS()
      tree.add_nodes_from(expSMS.nodes)
      tree.add_edges_from(expSMS.edgeIndices)
      tree.sort()

      nodeListB = []
      for node,nodeIndex in zip(tree.nodes,tree.nodeIndices):
          nodeListB.append((str(node),nodeIndex,
                          tree.nodeCanonName(nodeIndex)))
      self.assertEqual(nodesCanonNames,nodeListB)

      stringProc = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > b,W+,b), (H-(4) > W-,higgs), (H-(2) > W-,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      tree = TheorySMS()
      tree.add_nodes_from(expSMS.nodes)
      tree.add_edges_from(expSMS.edgeIndices)


      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('H0', 1), 
                                          ('H-', 2), ('H+', 3), 
                                          ('H-', 4), ('b', 5), 
                                          ('W+', 6), ('b', 7), 
                                          ('W-', 8), ('higgs', 9), 
                                          ('W-', 10), ('higgs', 11)])
      self.assertEqual(edges,[('H+', 'W+'), ('H+', 'b'), 
                              ('H+', 'b'), ('H-', 'W-'), 
                              ('H-', 'W-'), ('H-', 'higgs'), 
                              ('H-', 'higgs'), ('H0', 'H+'), 
                              ('H0', 'H-'), ('PV', 'H-'), 
                              ('PV', 'H0')])
      
      tree.sort()
      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('H-', 1), 
                                          ('H0', 2), ('higgs', 3), 
                                          ('W-', 4), ('H-', 5), 
                                          ('H+', 6), ('higgs', 7), 
                                          ('W-', 8), ('b', 9), 
                                          ('b', 10), ('W+', 11)])
      self.assertEqual(edges,[('H+', 'W+'), ('H+', 'b'), 
                              ('H+', 'b'), ('H-', 'W-'), 
                              ('H-', 'W-'), ('H-', 'higgs'), 
                              ('H-', 'higgs'), ('H0', 'H+'), 
                              ('H0', 'H-'), ('PV', 'H-'), 
                              ('PV', 'H0')])
      


      nodesCanonNames = [('PV', 0, 111010011101001101010000),
                        ('H-', 1, 110100),
                        ('H0', 2, 1110100110101000),
                        ('higgs', 3, 10),
                        ('W-', 4, 10),
                        ('H-', 5, 110100),
                        ('H+', 6, 11010100),
                        ('higgs', 7, 10),
                        ('W-', 8, 10),
                        ('b', 9, 10),
                        ('b', 10, 10),
                        ('W+', 11, 10)]
      nodeList = []
      for node,nodeIndex in zip(tree.nodes,tree.nodeIndices):
        nodeList.append((str(node),nodeIndex,tree.nodeCanonName(nodeIndex)))
      self.assertEqual(nodeList,nodesCanonNames)


      # Tree should be sorted after setting global properties:
      stringProc = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > W+,b,b), (H-(4) > W-,higgs), (H-(2) > W-,higgs)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      tree = TheorySMS()
      tree.add_nodes_from(expSMS.nodes)
      tree.add_edges_from(expSMS.edgeIndices)
      tree.maxWeight = 1.0*fb
      tree.prodXSec = 1.0*fb
      tree.setGlobalProperties()
      
      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('H-', 1), 
                                          ('H0', 2), ('higgs', 3), 
                                          ('W-', 4), ('H-', 5), 
                                          ('H+', 6), ('higgs', 7), 
                                          ('W-', 8), ('b', 9), 
                                          ('b', 10), ('W+', 11)])
      self.assertEqual(edges,[('H+', 'W+'), ('H+', 'b'), 
                              ('H+', 'b'), ('H-', 'W-'), 
                              ('H-', 'W-'), ('H-', 'higgs'), 
                              ('H-', 'higgs'), ('H0', 'H+'), 
                              ('H0', 'H-'), ('PV', 'H-'), 
                              ('PV', 'H0')])
      
      stringProc = "(PV > st_1~(1),st_1(2)), (st_1~(1) > N2~(3),t-), (st_1(2) > N2(4),t+), (N2~(3) > N1~,e-,e+), (N2(4) > N1,nu,nu)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      tree = TheorySMS()
      tree.add_nodes_from(expSMS.nodes)
      tree.add_edges_from(expSMS.edgeIndices)
      tree.maxWeight = 1.0*fb
      tree.prodXSec = 1.0*fb    
      tree.setGlobalProperties()

      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('st_1~', 1), 
                                          ('st_1', 2), ('t-', 3), 
                                          ('N2~', 4), ('t+', 5), 
                                          ('N2', 6), ('N1~', 7), 
                                          ('e-', 8), ('e+', 9), 
                                          ('N1', 10), ('nu', 11), 
                                          ('nu', 12)])
      self.assertEqual(edges,[('N2', 'N1'), ('N2', 'nu'), 
                              ('N2', 'nu'), ('N2~', 'N1~'), 
                              ('N2~', 'e+'), ('N2~', 'e-'), 
                              ('PV', 'st_1'), ('PV', 'st_1~'), 
                              ('st_1', 'N2'), ('st_1', 't+'), 
                              ('st_1~', 'N2~'), ('st_1~', 't-')])
      

      stringProc = "(PV > st_1~(1),st_1(2)), (st_1~(1) > N2~(3),t-), (st_1(2) > N2(4),t+), (N2~(3) > N1~,e+,e-), (N2(4) > N1,nu,nu)"
      # Hack to create a theory element from a string:
      expSMS = ExpSMS.from_string(stringProc, model=model)
      treeB = TheorySMS()
      treeB.add_nodes_from(expSMS.nodes)
      treeB.add_edges_from(expSMS.edgeIndices)
      treeB.maxWeight = 1.0*fb
      treeB.prodXSec = 1.0*fb

      nodes_and_indices = getNodesIndices(tree)
      edges = getEdges(tree)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('st_1~', 1), 
                                          ('st_1', 2), ('t-', 3), 
                                          ('N2~', 4), ('t+', 5), 
                                          ('N2', 6), ('N1~', 7), 
                                          ('e-', 8), ('e+', 9), 
                                          ('N1', 10), ('nu', 11), 
                                          ('nu', 12)])
      self.assertEqual(edges,[('N2', 'N1'), ('N2', 'nu'), 
                              ('N2', 'nu'), ('N2~', 'N1~'), 
                              ('N2~', 'e+'), ('N2~', 'e-'), 
                              ('PV', 'st_1'), ('PV', 'st_1~'), 
                              ('st_1', 'N2'), ('st_1', 't+'), 
                              ('st_1~', 'N2~'), ('st_1~', 't-')])
      
      self.assertTrue(tree == treeB)

if __name__ == "__main__":
    unittest.main()                         