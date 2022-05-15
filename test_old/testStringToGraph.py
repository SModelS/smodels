#!/usr/bin/env python3

"""
.. module:: testStringToGraph
   :synopsis: Tests the conversion of an element in bracket notation to a graph

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
import networkx as nx
from smodels.theory.auxiliaryFunctions import stringToGraph


class stringToGraphTest(unittest.TestCase):
    def testConversion(self):

        elementStr = "[[['e+','e-']],[['mu-'],['jet','jet'],['L']]]"
        gA = stringToGraph(elementStr)
        
        paths = [['PV','X','MET'], ['PV','X','e+'],['PV','X','e-'],
                 ['PV','Y','Z','W','MET'],['PV','Y','Z','W','L'],['PV','Y','Z','jetA'],['PV','Y','Z','jetB'],['PV','Y','mu-']]
        g,root = nx.prefix_tree(paths)
        g.remove_node(root)
        for n in g.nodes():
            if g.out_degree[n] == 0:
                g.remove_node(n)
                break
        
        #Set labels for comparison
        for n in g.nodes():
            if g.nodes[n]['source'] in ['X','Y','Z','W']:
                g.nodes[n]['label'] = 'anyOdd'
            elif g.nodes[n]['source'] in ['META','METB']:
                g.nodes[n]['label'] = 'MET'
            elif g.nodes[n]['source'] in ['jetA','jetB']:
                g.nodes[n]['label'] = 'jet'
            else:    
                g.nodes[n]['label'] = g.nodes[n]['source']        
        #Check if graphs have the same structure
        self.assertTrue(nx.isomorphism.is_isomorphic(g,gA)) 
        #Check if graphs have distinct labels:
        def compareNode(n1,n2):
            return  n1['label'] == n2['particle'].label
        self.assertTrue(nx.isomorphism.is_isomorphic(g,gA,node_match=compareNode))

if __name__ == "__main__":
    unittest.main()
