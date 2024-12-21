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
from smodels.base.smodelsLogging import setLogLevel
from smodels.base import runtime
from smodels.tools.particlesLoader import load
from smodels.base.smodelsLogging import logger
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList,proton,q
from smodels.experiment.defaultFinalStates import anyBSM, MET
from smodels.base.vertexGraph import VertexGraph
from smodels.base.inclusiveObjects import InclusiveValue

setLogLevel('error')


class VertexTest(unittest.TestCase):
    
    def testVertexCanonicalForm(self):

        filename = "./testFiles/slha/gluino_squarks.slha"
        runtime.modelFile = 'mssmQNumbers.slha'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)
        v = model.vertices[0]
        
        # Switching incoming and outgoing plus charge conjugation should
        # result in the same vertex
        v2 = VertexGraph(outgoing=[v.indexToNode(0).chargeConjugate()],
                 incoming=[p.chargeConjugate() for p in v.daughters(0)])
        self.assertEqual(v,v2)
        
        # Transposing one particle from the final state to the initial state
        # with charge conjugation should also give the same vertex
        v3 = VertexGraph(incoming=[v.indexToNode(0),v.daughters(0)[0].chargeConjugate()],
                 outgoing=[p for p in v.daughters(0)[1:]])
        self.assertEqual(v,v3)

        # Now charge-conjugation of only one state should not be equivalent
        v4 = VertexGraph(incoming=[v.indexToNode(0).chargeConjugate()],
                 outgoing=[p for p in v.daughters(0)])      
        self.assertNotEqual(v,v4)

    def testVertexConjugate(self):

        runtime.modelFile = 'mssmQNumbers.slha'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        go = model.pdgToParticle(1000021)
        x1m = model.pdgToParticle(-1000024)
        t = model.pdgToParticle(6)
        bbar = model.pdgToParticle(-5)
        v = VertexGraph(incoming=[go],outgoing=[x1m,t,bbar])
        # v = x1+ > go t+ b~
        x1p = model.pdgToParticle(1000024)
        tbar = model.pdgToParticle(-6)
        b = model.pdgToParticle(5)
        v2 = VertexGraph(incoming=[go],outgoing=[x1p,tbar,b])
        # v2 = x1+ > go t+ b~
        self.assertEqual(v,v2)
        v2Conj = v2.chargeConjugate()
        self.assertEqual(v,v2Conj)
        

    def testVertexMatching(self):

        filename = "./testFiles/slha/TRV1_1800_300_300.slha"
        runtime.modelFile = './testFiles/slha/TRV1_1800_300_300.slha'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)
        self.assertEqual(len(model.vertices),3)
        
        
        anyBSM.pdg = InclusiveValue() # Just needs to be larger than the proton pdg
        dm = model.pdgToParticle(52)
        dmbar = model.pdgToParticle(-52)
        y1 = model.pdgToParticle(55)
        prodV1 = VertexGraph(incoming=[proton,proton],outgoing=[anyBSM])
        nmatches = 0
        for v in model.vertices:
            if v.matchTo(prodV1) is not None:
                nmatches +=1
        self.assertEqual(nmatches,2)

        prodV2 = VertexGraph(incoming=[anyBSM],outgoing=[q,q])
        nmatches = 0
        for v in model.vertices:
            if v.matchTo(prodV2) is not None:
                nmatches +=1
        self.assertEqual(nmatches,2)

        prodV3 = VertexGraph(incoming=[y1],outgoing=[dm,dmbar])
        nmatches = 0
        for v in model.vertices:
            if v.matchTo(prodV3) is not None:
                nmatches +=1
        self.assertEqual(nmatches,1)

        prodV4 = VertexGraph(incoming=[dm,dmbar],outgoing=[y1.chargeConjugate()])
        nmatches = 0
        for v in model.vertices:
            if v.matchTo(prodV4) is not None:
                nmatches +=1
        self.assertEqual(nmatches,1)




if __name__ == "__main__":
    unittest.main()
