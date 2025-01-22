#!/usr/bin/env python3

"""
.. module:: testConvertNotation
   :synopsis: Testing Bracket to ProcessStr convertion

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.expAuxiliaryFuncs import bracketToProcessStr
from smodels.experiment.expSMS import ExpSMS
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.defaultFinalStates import finalStates


class TestConvertNotation(unittest.TestCase):
    
    # Tests that the function can handle a simple bracket notation
    def test_bracket_notation(self):
        stringEl = "[ [['e-','nu'], ['jet','jet'] ], [ ['L','nu'] ] ]"
        output = bracketToProcessStr(stringEl)

        self.assertEqual(output,'(PV(0) > anyBSM(1),anyBSM(2)), (anyBSM(1) > anyBSM(3),e-(4),nu(5)), (anyBSM(3) > MET(6),jet(7),jet(8)), (anyBSM(2) > MET(9),L(10),nu(11))')

        stringEl = "[ [['e-','nu'], ['jet','jet'] ], [ ['L','nu'] ] ]"
        output = bracketToProcessStr(stringEl,finalState=['MET','HSCP'],intermediateState=[['squark','gluino'],['squark']])
        self.assertEqual(output,'(PV(0) > squark(1),squark(2)), (squark(1) > gluino(3),e-(4),nu(5)), (gluino(3) > MET(6),jet(7),jet(8)), (squark(2) > HSCP(9),L(10),nu(11))')

        stringEl = "[ [ ['jet','jet'] ], [ ['L','nu'] ] ]"
        output = bracketToProcessStr(stringEl,finalState=['MET','HSCP'],intermediateState=[['gluino'],['anyBSM']])
        self.assertEqual(output,'(PV(0) > gluino(1),anyBSM(2)), (gluino(1) > MET(3),jet(4),jet(5)), (anyBSM(2) > HSCP(6),L(7),nu(8))')

    def test_convert_to_sms(self):

        procString = '(PV(0) > gluino(1),anyBSM(2)), (gluino(1) > MET(3),jet(4),jet(5)), (anyBSM(2) > HSCP(6),L(7),nu(8))'
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(procString, model=finalStates)
        tree = TheorySMS()
        tree.add_nodes_from(expSMS.nodes)
        tree.add_edges_from(expSMS.edgeIndices)

        PV = finalStates.getParticle(label='PV')
        gluino = finalStates.getParticle(label='gluino')
        anyBSM = finalStates.getParticle(label='anyBSM')
        MET = finalStates.getParticle(label='MET')
        jet = finalStates.getParticle(label='jet')
        HSCP = finalStates.getParticle(label='HSCP')
        L = finalStates.getParticle(label='L')
        nu = finalStates.getParticle(label='nu')
        self.assertEqual([node.particle for node in tree.nodes],[PV, gluino, anyBSM, MET, jet, jet, HSCP, L, nu])
        edges = [(mom.particle,daughter.particle) for mom,daughter in tree.edges]
        self.assertEqual(edges,[(PV, gluino), (PV, anyBSM), (gluino, MET), 
                                     (gluino, jet), (gluino, jet), (anyBSM, HSCP), 
                                     (anyBSM, L), (anyBSM, nu)])

    def test_convert_inclusive(self):      

        stringEl = "[[['*']],[]]"
        output = bracketToProcessStr(stringEl,finalState=['MET','MET'])
        self.assertEqual(output,'(PV(0) > anyBSM(1),MET(2)), (anyBSM(1) > MET(3),anySM(4))')

        stringEl = "[[['*','*']],[]]"
        output = bracketToProcessStr(stringEl,finalState=['MET','MET'])
        self.assertEqual(output,'(PV(0) > anyBSM(1),MET(2)), (anyBSM(1) > MET(3),anySM(4),anySM(5))')

        stringEl = "[['*'],[]]"
        output = bracketToProcessStr(stringEl,finalState=['MET','MET'])
        self.assertEqual(output,'(PV(0) > InclusiveNode(1),MET(2)), (InclusiveNode(1) > MET(3),*anySM(4))')

        stringEl = "[[['e-','e+']],['*']]"
        output = bracketToProcessStr(stringEl,finalState=['MET','HSCP']) 
        self.assertEqual(output,'(PV(0) > anyBSM(1),InclusiveNode(2)), (anyBSM(1) > MET(3),e+(4),e-(5)), (InclusiveNode(2) > HSCP(6),*anySM(7))')


if __name__ == "__main__":
    unittest.main()                         