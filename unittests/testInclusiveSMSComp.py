#!/usr/bin/env python3

"""
.. module:: testInclusiveSMSComp
   :synopsis: Testing comparison of inclusive SMS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.experiment.defaultFinalStates import finalStates
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV
from unitTestHelpers import getNodesIndices, getEdges

slhafile="./testFiles/slha/lightEWinos.slha"
model = Model(BSMList,SMList)
model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)


class TestInclusiveSMSComp(unittest.TestCase):

    
    
    def test_comp(self):

        smsA = ExpSMS.from_string("[ [ ['e-','nu'] ], [['ta+','ta-'],['u,u~']] ]",model=model,
              intermediateState=[['C1-'],['N2','gluino']], finalState=['N1','N1'])
        smsB = ExpSMS.from_string("[ ['*'], [ ['nu','L'] ] ]",model=finalStates,finalState=['MET','MET'])

        matched = smsA.matchesTo(smsB)
        nodes_and_indices = getNodesIndices(matched)
        edges = getEdges(matched)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('anyBSM', 1), 
                                            ('Inclusive', 2), ('MET', 3), 
                                            ('L', 4), ('nu', 5), 
                                            ('MET', 9), ('*anySM', 11)])
        self.assertEqual(edges,[('Inclusive', '*anySM'), 
                                ('Inclusive', 'MET'), 
                                ('PV', 'Inclusive'), ('PV', 'anyBSM'), 
                                ('anyBSM', 'L'), ('anyBSM', 'MET'), 
                                ('anyBSM', 'nu')])
        
        matched = smsB.matchesTo(smsA)
        nodes_and_indices = getNodesIndices(matched)
        edges = getEdges(matched)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('N2', 1), 
                                            ('C1-', 2), ('N1', 3), 
                                            ('q', 4), ('N1', 5), 
                                            ('e-', 6), ('nu', 7), 
                                            ('ta-', 8), ('q', 9), 
                                            ('gluino', 10), ('ta+', 11)])
        
        self.assertEqual(edges,[('C1-', 'N1'), ('C1-', 'e-'), 
                                ('C1-', 'nu'), ('N2', 'gluino'), 
                                ('N2', 'ta+'), ('N2', 'ta-'), 
                                ('PV', 'C1-'), ('PV', 'N2'), 
                                ('gluino', 'N1'), ('gluino', 'q'), 
                                ('gluino', 'q')])


        smsC = ExpSMS.from_string("[ [ ['e-','q'] ], [['e+','ta-']] ]",model=model,
              intermediateState=[['C1-'],['N2']], finalState=['N1','C2+'])

        matched = smsB.matchesTo(smsC)
        self.assertIs(matched,None)
        matched = smsC.matchesTo(smsB)
        self.assertIs(matched,None)
        matched = smsC.matchesTo(smsA)
        self.assertIs(matched,None)
        matched = smsA.matchesTo(smsC)
        self.assertIs(matched,None)


        smsA = ExpSMS.from_string("[ [ ['*','*'] ], ['*'] ]",
                                  model=finalStates, finalState=['HSCP','MET'])
        smsB = ExpSMS.from_string("[ ['*'], [ ['*'] ] ]",
                                  model=finalStates, finalState=['MET','HSCP'])
        
        matched = smsA.matchesTo(smsB)
        self.assertIs(matched,None)
        matched = smsB.matchesTo(smsA)
        self.assertIs(matched,None)


        smsA = ExpSMS.from_string("[ [ ['*','*'] ], ['*'] ]",model=finalStates, 
              finalState=['HSCP','MET'])
        smsB = ExpSMS.from_string("[ [ ['*'] ], ['*'] ]",model=finalStates,
              finalState=['MET','HSCP'])
        
        matched = smsA.matchesTo(smsB)
        nodes_and_indices = getNodesIndices(matched)
        edges = getEdges(matched)

        self.assertEqual(nodes_and_indices,[('PV', 0), ('Inclusive', 1),
                                             ('anyBSM', 2), ('HSCP', 3),
                                            ('*anySM', 5), ('MET', 6), 
                                            ('anySM', 7)])
        
        self.assertEqual(edges,[('Inclusive', '*anySM'),
                                 ('Inclusive', 'HSCP'),
                                 ('PV', 'Inclusive'), ('PV', 'anyBSM'), 
                                 ('anyBSM', 'MET'), ('anyBSM', 'anySM')])
        

        matched = smsB.matchesTo(smsA)
        nodes_and_indices = getNodesIndices(matched)
        edges = getEdges(matched)
        
        self.assertEqual(nodes_and_indices,[('PV', 0), ('Inclusive', 1),
                                             ('anyBSM', 2), ('MET', 3),
                                            ('*anySM', 4), ('HSCP', 5), 
                                            ('anySM', 6), ('anySM', 7)])
        
        self.assertEqual(edges,[('Inclusive', '*anySM'), 
                                ('Inclusive', 'MET'), ('PV', 'Inclusive'), 
                                ('PV', 'anyBSM'), ('anyBSM', 'HSCP'), 
                                ('anyBSM', 'anySM'), ('anyBSM', 'anySM')])

            
if __name__ == "__main__":
    unittest.main()                         