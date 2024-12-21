#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
import itertools
from smodels.base.smodelsLogging import setLogLevel
from smodels.base import runtime
from smodels.tools.particlesLoader import load
from smodels.base.smodelsLogging import logger
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.base.productionGraph import get2to1ProdGraph,get2to2ProdGraph
from smodels.base.vertexGraph import VertexGraph


setLogLevel('error')


class ProdGraphTest(unittest.TestCase):
    
    def test2to1(self):
        filename = "./testFiles/slha/TRV1_1800_300_300.slha"
        runtime.modelFile = './testFiles/slha/TRV1_1800_300_300.slha'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)
        matches = get2to1ProdGraph(model)

        self.assertEqual(len(matches),2)

        u = model.pdgToParticle(2)
        d = model.pdgToParticle(1)
        y1 = model.pdgToParticle(55)
        v1 = VertexGraph(incoming=[u,u.chargeConjugate()], outgoing=[y1])
        v2 = VertexGraph(incoming=[d,d.chargeConjugate()],outgoing=[y1])
        
        self.assertTrue(v1 in matches)
        self.assertTrue(v2 in matches)

        # Now require only the u u~ initial state
        matches = get2to1ProdGraph(model,initialStates=[u,u.chargeConjugate()],finalState=y1)
        self.assertEqual(len(matches),1)
        self.assertEqual(matches[0],v1)

        # Finally ask for a process which is not possible
        matches = get2to1ProdGraph(model,initialStates=[u,u],finalState=y1)
        self.assertEqual(len(matches),0)

    def test2to2(self):
        filename = "./testFiles/slha/TRV1_1800_300_300.slha"
        runtime.modelFile = './testFiles/slha/TRV1_1800_300_300.slha'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)
        u = model.pdgToParticle(2)
        ubar = model.pdgToParticle(-2)
        d = model.pdgToParticle(1)
        dbar = model.pdgToParticle(-1)
        y1 = model.pdgToParticle(55)
        xd = model.pdgToParticle(52)
        xdbar = model.pdgToParticle(-52)

        matchesDict = get2to2ProdGraph(model)
        self.assertEqual(sorted(matchesDict.keys()),sorted(['s','t', 'u']))

        self.assertEqual(len(matchesDict['s']),2)
        self.assertEqual(len(matchesDict['t']),2)
        self.assertEqual(len(matchesDict['u']),2)

        #s-channel:
        vR = VertexGraph(incoming=[y1],outgoing=[xd,xdbar])
        v1L = VertexGraph(incoming=[u,ubar],outgoing=[y1])
        v2L = VertexGraph(incoming=[d,dbar],outgoing=[y1])
        for va,vb in itertools.product([v1L,v2L],[vR]):
            self.assertTrue((va,vb) in matchesDict['s'])

        #t-channel and u-channel:
        v1L = VertexGraph(incoming=[y1,ubar],outgoing=[ubar])
        v1R = VertexGraph(incoming=[ubar],outgoing=[y1,ubar])
        v2L = VertexGraph(incoming=[y1,dbar],outgoing=[dbar])
        v2R = VertexGraph(incoming=[dbar],outgoing=[y1,dbar])
        for va,vb in zip([v1L,v2L],[v1R,v2R]):
            self.assertTrue((va,vb) in matchesDict['t'])
            self.assertTrue((va,vb) in matchesDict['u'])


if __name__ == "__main__":
    unittest.main()
