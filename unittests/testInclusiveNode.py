#!/usr/bin/env python3

"""
.. module:: testInclusiveNodeClass
   :synopsis: Tests the theory.particleNode.InclusiveParticleNode class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,time
sys.path.append('../')
from smodels.theory.particleNode import InclusiveParticleNode
from smodels.experiment.expAuxiliaryFuncs import bracketToProcessStr
from smodels.experiment.defaultFinalStates import finalStates
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.tools.inclusiveObjects import InclusiveValue
from unitTestHelpers import theorySMSFromString as fromString
import networkx as nx
import unittest


class InclusiveNodeTest(unittest.TestCase):

    def testString(self):

        stringEl = "[ [ ['L','nu'] ], ['*'] ]"
        output = bracketToProcessStr(stringEl)
        self.assertEqual(output.replace(" ",""),
                         "(PV>anyBSM(1),InclusiveNode(2)),(anyBSM(1)>L,nu,MET),(InclusiveNode(2)>anySM,MET)")
        stringEl = "[ [['*'],['e+']], [ ['L','nu'] ] ]"
        output2 = bracketToProcessStr(stringEl)
        self.assertEqual(output2,"(PV > anyBSM(1),anyBSM(3)),(anyBSM(1) > anySM,anyBSM(2)),(anyBSM(2) > e+,MET),(anyBSM(3) > L,nu,MET)")

        stringEl = "[ ['*'], [ ['L','nu'] ] ]"
        output = bracketToProcessStr(stringEl)
        self.assertEqual(output,"(PV > InclusiveNode(1),anyBSM(2)),(InclusiveNode(1) > anySM,MET),(anyBSM(2) > L,nu,MET)")

    def testGraph(self):

        stringEl = "[ [ ['L','nu'] ], ['*'] ]"
        T = ExpSMS.from_string(stringEl,model=finalStates)
        nodes = [str(n) for n in T.nodes]
        self.assertEqual(nodes,['PV', 'anyBSM', 'Inclusive', 'L', 'nu', 'MET'])
        edges = [(str(edge[0]),str(edge[1])) for edge in T.edges]
        self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'Inclusive'),
                                ('anyBSM', 'L'), ('anyBSM', 'nu'), ('anyBSM', 'MET')])

        stringEl = "[ [['*'],['e+']], [ ['L','nu'] ] ]"
        output2 = bracketToProcessStr(stringEl)
        procString = output2
        T = ExpSMS.from_string(procString,model=finalStates)
        nodes = [str(n) for n in T.nodes]
        self.assertEqual(nodes,['PV', 'anyBSM', 'anyBSM',
                                'anyBSM', 'anySM', 'e+', 'MET', 'L', 'nu', 'MET'])
        edgeIndices = [(str(edge[0]),str(edge[1])) for edge in T.edges]
        self.assertEqual(T.edgeIndices,[(0, 1), (0, 3), (1, 2),
                         (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (3, 9)])


    def testCanonName(self):

        stringEl = "[ [ ['L','nu'] ], ['*'] ]"
        T = ExpSMS.from_string(stringEl,model=finalStates)
        cNamesDict = {'PV' : InclusiveValue(),'anyBSM' : 11010100, 'Inclusive': InclusiveValue(),
                        'L' : 10, 'nu' : 10, 'MET': 10}
        for nodeIndex in T.nodeIndices:
            node = T.indexToNode(nodeIndex)
            cname = T.nodeCanonName(nodeIndex)
            self.assertEqual(cname,cNamesDict[str(node)])
            self.assertEqual(str(cname),str(cNamesDict[str(node)]))

    def testComparison(self):


        slhafile="../inputFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV)

        treeA = fromString("[ [ ['e-','nu'] ], [['ta+','ta-'],['u,u~']] ]",model=model,
                      intermediateState=[['C1-'],['N2','gluino']], finalState=['N1','N1'])

        # Experiment SMS
        treeB = ExpSMS.from_string("[ ['*'], [ ['L','nu'] ] ]",model=finalStates)

        treeC = fromString("[ [ ['e-','q'] ], [['e+','ta-']] ]",model=model,
                            intermediateState=[['C1-'],['N2','gluino']],
                            finalState=['N1','C2+'])

        matchedEl = treeB.matchesTo(treeA)
        nodes = [str(n) for n in matchedEl.nodes]
        self.assertEqual(nodes,['PV', 'Inclusive', 'C1-', 'e-', 'nu', 'N1'])
        edges = [(str(edge[0]),str(edge[1])) for edge in matchedEl.edges]
        self.assertEqual(edges,[('PV', 'Inclusive'), ('PV', 'C1-'),
                                ('C1-', 'e-'), ('C1-', 'nu'), ('C1-', 'N1')])

        self.assertTrue(treeB.matchesTo(treeC) is None)

        mass = [None, None, 1.34E+02, 5.00E-04, 0.00E+00, 6.81E+01]
        for im,m in enumerate(matchedEl.mass):
            if m is None:
                self.assertEqual(m,mass[im])
            else:
                self.assertEqual(float('%1.2e' %(m.asNumber(GeV))),mass[im])


        stringEl = "[ [['*','t-']], [ ['L','nu'] ] ]"
        sms = ExpSMS.from_string(stringEl,model=finalStates)

        stringEl = "[ [['t-','nu']], [ ['e-','nu'] ] ]"
        smsB = ExpSMS.from_string(stringEl,model=model,finalState=['N1','N1'],intermediateState=[['gluino'],['gluino']])

        matchedEl = sms.matchesTo(smsB)
        nodes = [str(n) for n in matchedEl.nodes]
        self.assertEqual(nodes,['PV', 'gluino', 'gluino', 'nu', 't-', 'N1', 'e-', 'nu', 'N1'])


        stringEl = "[ [['t-','nu']], [ ['e-','e+'] ] ]"
        smsC = ExpSMS.from_string(stringEl,model=model,finalState=['N1','N1'],intermediateState=[['gluino'],['gluino']])
        self.assertTrue(sms.matchesTo(smsC) is None)




if __name__ == "__main__":
    unittest.main()
