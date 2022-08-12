#!/usr/bin/env python3

"""
.. module:: testInclusiveNodeClass
   :synopsis: Tests the theory.particleNode.InclusiveParticleNode class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,time
sys.path.append('../')
from smodels.base.particleNode import InclusiveParticleNode
from smodels.experiment.expAuxiliaryFuncs import bracketToProcessStr
from smodels.experiment.defaultFinalStates import finalStates
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.base.inclusiveObjects import InclusiveValue
from unitTestHelpers import theorySMSFromString as fromString
import networkx as nx
import unittest


class InclusiveNodeTest(unittest.TestCase):

    def testString(self):

        stringEl = "[ [ ['L','nu'] ], ['*'] ]"
        output = bracketToProcessStr(stringEl)
        self.assertEqual(output.replace(" ",""),
                         "(PV>anyBSM(1),InclusiveNode(2)),(anyBSM(1)>L,nu,MET),(InclusiveNode(2)>*anySM,MET)")
        stringEl = "[ [['*'],['e+']], [ ['L','nu'] ] ]"
        output2 = bracketToProcessStr(stringEl)
        self.assertEqual(output2,"(PV > anyBSM(1),anyBSM(3)),(anyBSM(1) > anySM,anyBSM(2)),(anyBSM(2) > e+,MET),(anyBSM(3) > L,nu,MET)")

        stringEl = "[ ['*'], [ ['L','nu'] ] ]"
        output = bracketToProcessStr(stringEl)
        self.assertEqual(output,"(PV > InclusiveNode(1),anyBSM(2)),(InclusiveNode(1) > *anySM,MET),(anyBSM(2) > L,nu,MET)")

    def testGraph(self):

        stringEl = "[ [ ['L','nu'] ], ['*'] ]"
        T = ExpSMS.from_string(stringEl,model=finalStates)
        nodes = [str(n) for n in T.nodes]
        self.assertEqual(nodes,['PV', 'anyBSM', 'Inclusive', 'L', 'nu', 'MET','*anySM','MET'])
        edges = sorted([(str(edge[0]),str(edge[1])) for edge in T.edges])
        self.assertEqual(edges,sorted([('PV', 'anyBSM'), ('PV', 'Inclusive'),
                                       ('Inclusive','*anySM'), ('Inclusive','MET'),
                                       ('anyBSM', 'L'), ('anyBSM', 'nu'), ('anyBSM', 'MET')]))

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
                        'L' : 10, 'nu' : 10, 'MET': 10, '*anySM' : InclusiveValue()}
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
        for nodeIndex,node in zip(treeB.nodeIndices,treeB.nodes):
            self.assertEqual(node,matchedEl.indexToNode(nodeIndex))

        nodesDict = {0 : 'PV', 1 : 'N2', 2 : 'C1-', 3 : 'ta-', 4 : 'N1', 5 : 'e-', 6 : 'nu', 7 : 'N1'}
        for nodeIndex,nodeStr in nodesDict.items():
            self.assertEqual(str(matchedEl.indexToNode(nodeIndex)),nodeStr)

        nodes = [str(n) for n in matchedEl.nodes]
        self.assertEqual(nodes,['PV', 'N2', 'C1-', 'ta-', 'N1', 'e-', 'nu', 'N1', 'gluino', 'q', 'q', 'ta+'])
        edges = [(str(edge[0]),str(edge[1])) for edge in matchedEl.edges]
        self.assertEqual(edges,[('PV', 'N2'), ('PV', 'C1-'), ('N2','ta-'),
                                ('N2','gluino'), ('N2','ta+'), ('C1-', 'e-'), ('C1-', 'nu'), ('C1-', 'N1'),
                                ('gluino','N1'), ('gluino','q'),('gluino','q')])

        self.assertTrue(treeB.matchesTo(treeC) is None)

        n2 = model.getParticle(label='N2')
        n1 = model.getParticle(label='N1')
        c1 = model.getParticle(label='C1-')
        mass = {1 : n2.mass, 2 : c1.mass, 4 : n1.mass, 7 : n1.mass}
        for nodeIndex in treeB.nodeIndices:
            node = matchedEl.indexToNode(nodeIndex)
            if node.isSM:
                continue
            m = float('%1.2e' %(node.mass.asNumber(GeV)))
            m_default = float('%1.2e' %(mass[nodeIndex].asNumber(GeV)))
            self.assertEqual(m,m_default)


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
