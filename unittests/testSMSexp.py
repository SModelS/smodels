#!/usr/bin/env python3

"""
.. module:: testSMSexp
   :synopsis: Testing creation of ExpSMS objects

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
from smodels.experiment.defaultFinalStates import finalStates
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV
from unitTestHelpers import getNodesIndices, getEdges



class TestSMSExp(unittest.TestCase):
    
    
    def test_create_from_bracket(self):

      sms = "[[['e+','e-'],['jet','jet']],[['mu-','nu'],['ta+']]]"
      sms = ExpSMS.from_string(sms,model=finalStates)
      sms.bfs_sort(numberNodes=True)
      self.assertEqual(sms.canonName,1110101101000110101101010000)
      self.assertEqual(str(sms),'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),e+,e-), (anyBSM(2) > anyBSM(6),mu-,nu), (anyBSM(3) > MET,jet,jet), (anyBSM(6) > MET,ta+)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('anyBSM', 3), 
                                          ('e+', 4), ('e-', 5), 
                                          ('anyBSM', 6), ('mu-', 7), 
                                          ('nu', 8), ('MET', 9), 
                                          ('jet', 10), ('jet', 11), 
                                          ('MET', 12), ('ta+', 13)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'MET'), 
                              ('anyBSM', 'anyBSM'), ('anyBSM', 'anyBSM'), 
                              ('anyBSM', 'e+'), ('anyBSM', 'e-'), 
                              ('anyBSM', 'jet'), ('anyBSM', 'jet'), 
                              ('anyBSM', 'mu-'), ('anyBSM', 'nu'), 
                              ('anyBSM', 'ta+')])
      

      sms = "[ [['e-','nu'], ['jet','jet'] ], [ ['L','nu'] ] ]"
      fs=['MET','HSCP']
      intermediateStates=[['anyBSM','gluino'],['anyBSM']]
      sms = ExpSMS.from_string(sms,finalState=fs,
                              intermediateState=intermediateStates, 
                              model=finalStates)
      sms.bfs_sort(numberNodes=True)
      self.assertEqual(sms.canonName,111010100110101101010000)
      self.assertEqual(str(sms),'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > gluino(3),e-,nu), (anyBSM(2) > HSCP,L,nu), (gluino(3) > MET,jet,jet)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('gluino', 3), 
                                          ('e-', 4), ('nu', 5), 
                                          ('HSCP', 6), ('L', 7), 
                                          ('nu', 8), ('MET', 9), 
                                          ('jet', 10), ('jet', 11)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'HSCP'), ('anyBSM', 'L'), 
                              ('anyBSM', 'e-'), ('anyBSM', 'gluino'), 
                              ('anyBSM', 'nu'), ('anyBSM', 'nu'), 
                              ('gluino', 'MET'), ('gluino', 'jet'), 
                              ('gluino', 'jet')])


    def test_create_from_str(self):


      strProc = '(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > gluino(3),e-,nu), (anyBSM(2) > HSCP,L,nu), (gluino(3) > MET,jet,jet)'
      sms = ExpSMS.from_string(strProc, model=finalStates)

      self.assertEqual(sms.canonName,111010100110101101010000)
      self.assertEqual(str(sms),'(PV > anyBSM(1),anyBSM(2)), (anyBSM(1) > gluino(3),e-,nu), (anyBSM(2) > HSCP,L,nu), (gluino(3) > MET,jet,jet)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('gluino', 3), 
                                          ('e-', 4), ('nu', 5), 
                                          ('HSCP', 6), ('L', 7), 
                                          ('nu', 8), ('MET', 9), 
                                          ('jet', 10), ('jet', 11)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), 
                              ('anyBSM', 'HSCP'), ('anyBSM', 'L'), 
                              ('anyBSM', 'e-'), ('anyBSM', 'gluino'), 
                              ('anyBSM', 'nu'), ('anyBSM', 'nu'), 
                              ('gluino', 'MET'), ('gluino', 'jet'), 
                              ('gluino', 'jet')])
      
      strProc = "(PV > anyBSM(1)), (anyBSM(1) > e+, e-)"
      sms = ExpSMS.from_string(strProc, model=finalStates)
      self.assertEqual(sms.canonName,11101000)
      self.assertEqual(str(sms),'(PV > anyBSM(1)), (anyBSM(1) > e+,e-)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('e+', 2), ('e-', 3)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('anyBSM', 'e+'), 
                              ('anyBSM', 'e-')])
      

      strProc = "(PV > anyBSM(1)), (anyBSM(1) > anyBSM(2), higgs), (anyBSM(2) > jet,jet)"
      sms = ExpSMS.from_string(strProc, model=finalStates)
      self.assertEqual(sms.canonName,111011010000)
      self.assertEqual(str(sms),'(PV > anyBSM(1)), (anyBSM(1) > anyBSM(2),higgs), (anyBSM(2) > jet,jet)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('anyBSM', 2), ('higgs', 3), 
                                          ('jet', 4), ('jet', 5)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('anyBSM', 'anyBSM'), 
                              ('anyBSM', 'higgs'), ('anyBSM', 'jet'), 
                              ('anyBSM', 'jet')])
      
      strProc = "(PV > q,anyBSM(1)), (anyBSM(1) > e+, jet)"
      sms = ExpSMS.from_string(strProc,model=finalStates)
      self.assertEqual(sms.canonName,1101101000)
      self.assertEqual(str(sms),'(PV > anyBSM(1),q), (anyBSM(1) > e+,jet)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                          ('q', 2), ('e+', 3), 
                                          ('jet', 4)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'q'), 
                              ('anyBSM', 'e+'), ('anyBSM', 'jet')])
      

      strProc = '(PV > gluino(1),anyBSM(2)), (gluino(1) > anyBSM(3),anyBSM(4),jet), (anyBSM(2) > HSCP,L,nu), (anyBSM(3) > nu,nu), (anyBSM(4) > MET,e-)'
      sms = ExpSMS.from_string(strProc, model=finalStates)
      self.assertEqual(sms.canonName,11101010011011010011010000)
      self.assertEqual(str(sms),'(PV > gluino(1),anyBSM(2)), (gluino(1) > anyBSM(3),anyBSM(4),jet), (anyBSM(2) > HSCP,L,nu), (anyBSM(3) > nu,nu), (anyBSM(4) > MET,e-)')

      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('gluino', 1), 
                                          ('anyBSM', 2), ('anyBSM', 3), 
                                          ('anyBSM', 4), ('jet', 5), 
                                          ('HSCP', 6), ('L', 7), 
                                          ('nu', 8), ('nu', 9), 
                                          ('nu', 10), ('MET', 11), 
                                          ('e-', 12)])
      self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'gluino'), 
                              ('anyBSM', 'HSCP'), ('anyBSM', 'L'), 
                              ('anyBSM', 'MET'), ('anyBSM', 'e-'), 
                              ('anyBSM', 'nu'), ('anyBSM', 'nu'), 
                              ('anyBSM', 'nu'), ('gluino', 'anyBSM'), 
                              ('gluino', 'anyBSM'), ('gluino', 'jet')])
      

    def test_create_using_model(self):
       
      slhafile = './testFiles/slha/lightEWinos.slha'
      model = Model(BSMparticles=BSMList, SMparticles=SMList)
      model.updateParticles(inputFile=slhafile)

      strProc = "(PV > su_L(1),su_L~(2)), (su_L(1) > gluino(3),e-,nu), (su_L~(2) > C1+,d), (gluino(3) > N1,u,u~)"
      sms = ExpSMS.from_string(strProc, model=model)
      nodes_and_indices = getNodesIndices(sms)
      edges = getEdges(sms)

      self.assertEqual(nodes_and_indices,[('PV', 0), ('su_L', 1), 
                                          ('su_L~', 2), ('gluino', 3), 
                                          ('e-', 4), ('nu', 5), 
                                          ('C1+', 6), ('q', 7), 
                                          ('N1', 8), ('q', 9), 
                                          ('q', 10)])
      self.assertEqual(edges,[('PV', 'su_L'), ('PV', 'su_L~'), 
                              ('gluino', 'N1'), ('gluino', 'q'), 
                              ('gluino', 'q'), ('su_L', 'e-'), 
                              ('su_L', 'gluino'), ('su_L', 'nu'), 
                              ('su_L~', 'C1+'), ('su_L~', 'q')])
      

    def test_convert_bracket(self):
       
      slhafile = './testFiles/slha/lightEWinos.slha'
      model = Model(BSMparticles=BSMList, SMparticles=SMList)
      model.updateParticles(inputFile=slhafile)

      strProc = "(PV > su_L(1),su_L~(2)), (su_L(1) > gluino(3),e-,nu), (su_L~(2) > C1+,d), (gluino(3) > N1,u,u~)"
      sms = ExpSMS.from_string(strProc, model=model)
      smsB,finalState,intermediateState = sms.treeToBrackets()
      self.assertEqual(smsB,[[['e-', 'nu'], ['q', 'q']], [['q']]])
      self.assertEqual(finalState,['N1', 'C1+'])
      self.assertEqual(intermediateState,[['su_L', 'gluino'], ['su_L~']])


      strProc = "(PV > N1,C1+(1)), (C1+(1) > N1,u,d)"
      sms = ExpSMS.from_string(strProc,model=model)
      smsB,finalState,intermediateState = sms.treeToBrackets()
      self.assertEqual(smsB,[[['q', 'q']], []])
      self.assertEqual(finalState,['N1', 'N1'])
      self.assertEqual(intermediateState,[['C1+'], []])


    def test_getattr(self):
       
      slhafile = './testFiles/slha/lightEWinos.slha'
      model = Model(BSMparticles=BSMList, SMparticles=SMList)
      model.updateParticles(inputFile=slhafile)
      strProc = "(PV > N1,C1+(1)), (C1+(1) > N1,u,d)"
      sms = ExpSMS.from_string(strProc,model=model)

      self.assertEqual(sms.getFinalStates(),[3, 4, 5, 2])
      fs = [str(n) for n in sms.indexToNode(sms.getFinalStates())]
      self.assertEqual(fs,['N1', 'q', 'q', 'N1'])


if __name__ == "__main__":
    unittest.main()                         