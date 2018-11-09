#!/usr/bin/env python3

"""
.. module:: testStringToGraph
   :synopsis: Tests the conversion of an element in bracket notation to a graph

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from graph_tool import Graph, Vertex, topology
from smodels.experiment.finalStateParticles import finalStates
from smodels.theory.particle import Particle
from smodels.theory.auxiliaryFunctions import stringToGraph


class stringToGraphTest(unittest.TestCase):
    def testConversion(self):

        elementStr = "[[['e+','e-']],[['mu-'],['jet','jet'],['L']]]"
        gA = stringToGraph(elementStr)
        
        ep = finalStates.getParticlesWith(label='e+')[0]
        em = finalStates.getParticlesWith(label='e-')[0]
        mu = finalStates.getParticlesWith(label='mu-')[0]
        jet = finalStates.getParticlesWith(label='jet')[0]
        L = finalStates.getParticlesWith(label='L')[0]
        MET = finalStates.getParticlesWith(label='MET')[0]
        odd = finalStates.getParticlesWith(label='anyOdd')[0]
        
        
        # v0 -> v1 + v2, v1 -> (e-,e+,MET), v2 -> (v3,mu-), v3->(v4,jet,jet), v4->(L,MET) 
        g = Graph(directed=True)
        g.vertex_properties['particle'] = g.new_vertex_property("object")
        g.vertex_properties['label'] = g.new_vertex_property("string")
        g.vertex_properties['inclusive'] = g.new_vertex_property("bool",False)
        v0 = g.add_vertex()
        g.vp.particle[v0] = Particle(label='PV')
        v1 = g.add_vertex()
        g.vp.particle[v1] = odd
        g.add_edge(v0,v1)
        emV,epV,metV = g.add_vertex(3)
        g.vp.particle[emV] = em
        g.vp.particle[epV] = ep
        g.vp.particle[metV] = MET
        g.add_edge(v1,emV)
        g.add_edge(v1,epV)
        g.add_edge(v1,metV)
        
        
        v2 = g.add_vertex()
        g.vp.particle[v2] = odd
        g.add_edge(v0,v2)
        v3,muV = g.add_vertex(2)
        g.vp.particle[v3] = odd
        g.vp.particle[muV] = mu
        g.add_edge(v2,v3)
        g.add_edge(v2,muV)
        v4,jVa,jVb = g.add_vertex(3)
        g.vp.particle[v4] = odd
        g.vp.particle[jVa] = jet
        g.vp.particle[jVb] = jet
        g.add_edge(v3,v4)
        g.add_edge(v3,jVa)
        g.add_edge(v3,jVb)
        LV,metV = g.add_vertex(2)
        g.vp.particle[LV] = L
        g.vp.particle[metV] = MET
        g.add_edge(v4,LV)
        g.add_edge(v4,metV)
        
        #Set object IDs for comparison:
        for iv in range(g.num_vertices()):
            g.vp.label[iv] = g.vp.particle[iv].label
        
        self.assertTrue(topology.isomorphism(g,gA)) #Check if graphs have the same structure
        isomaps = topology.subgraph_isomorphism(g,gA,vertex_label=(g.vp.label,gA.vp.label))
        self.assertTrue(len(isomaps) == 2) #Check if graphs have the same labels (two mappings are possible corresponding to jet <-> jet)

if __name__ == "__main__":
    unittest.main()
