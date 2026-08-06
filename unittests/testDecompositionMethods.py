#!/usr/bin/env python3

"""
.. module:: testDecompositionMethods
   :synopsis: Testing methods used for decomposition

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
import numpy as np
from smodels.decomposition.theorySMS import TheorySMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition.decomposer import get_particle_order_dict, get_lightweight_decays, build_subtree_cacheFor, decayTupleObj
from smodels.base.particleNode import ParticleNode
from smodels.base.physicsUnits import fb, GeV


slhafile = './testFiles/slha/lightEWinos.slha'
model = Model(BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile)

pv = ParticleNode(model.getParticlesWith(label='PV')[0])
n1 = ParticleNode(model.getParticlesWith(pdg=1000022)[0])
n2 = ParticleNode(model.getParticlesWith(pdg=1000023)[0])
c2 = ParticleNode(model.getParticlesWith(pdg=1000037)[0])

def getNdecays(model,pdg):
    ptc = model.getParticlesWith(pdg=pdg)
    nDecays = 0
    if not ptc:
        return 0
    ptc = ptc[0]
    if not hasattr(ptc,'decays'):
        return 1
    if not ptc.decays:
        return 1
    if ptc.isStable():
        return 1
    
    for dec in ptc.decays:
        daughterPDG = [abs(pid) for pid in dec.ids if abs(pid) > 10000][0]
        nDecays += getNdecays(model,daughterPDG)
    
    return nDecays


class TestDecompositionMethods(unittest.TestCase):

    
    # Tests that the function can handle a simple bracket notation
    def test_decay_trees(self):
        particleOrderDict = get_particle_order_dict(model)
        decaysDict = get_lightweight_decays(model, particleOrderDict)

        n2_id = hash(n2.particle)
        decs = decaysDict[n2_id]

        labelDict = {hash(p): str(p) for p in model.SMparticles + model.BSMparticles}
        default_output = [['N2', 'N1', 'q', 'q',  0.187809422 ],
                    ['N2', 'N1', 'q', 'q',  0.187809422 ],
                    ['N2', 'N1', 'b', 'b',  0.180233089 ],
                    ['N2', 'N1', 'q', 'q',  0.144875155 ],
                    ['N2', 'N1', 'c', 'c',  0.144875155 ],
                    ['N2', 'N1', 'nu', 'nu',  0.0287784537 ],
                    ['N2', 'N1', 'nu', 'nu',  0.0287784537 ],
                    ['N2', 'N1', 'nu', 'nu',  0.0285446189 ],
                    ['N2', 'N1', 'ta+', 'ta-',  0.0219195363 ],
                    ['N2', 'N1', 'e+', 'e-',  0.0210071956 ],
                    ['N2', 'N1', 'mu+', 'mu-',  0.0210071956 ],
                    ['N2', 'N1', 'photon',  0.00436171735 ],
                    ['N2', 'C1+', 'q', 'q',  1.077496e-07 ],
                    ['N2', 'C1-', 'q', 'q',  1.077496e-07 ],
                    ['N2', 'C1+', 'c', 'q',  1.077496e-07 ],
                    ['N2', 'C1-', 'c', 'q',  1.077496e-07 ],
                    ['N2', 'C1+', 'e-', 'nu',  3.84522575e-08 ],
                    ['N2', 'C1-', 'e+', 'nu',  3.84522575e-08 ],
                    ['N2', 'C1+', 'mu-', 'nu',  3.84522575e-08 ],
                    ['N2', 'C1-', 'mu+', 'nu',  3.84522575e-08 ]]
        default_output = sorted(default_output)
        # Convert lightweight decay tuples (hash IDs) to labels for comparison:
        decs = sorted([[labelDict[dec.mom]]
                       + sorted([labelDict[d] for d in dec.daughters])
                       + [dec.br] for dec in decs])
        self.assertEqual(len(default_output),len(decs))
        for idec,dec in enumerate(decs):
            default_dec = default_output[idec]
            self.assertEqual(dec[:-1],default_dec[:-1])
            self.assertTrue(np.isclose(dec[-1],default_dec[-1],atol=0,rtol=1e-4))



    def test_sms_decomp(self):

        tree = TheorySMS()
        tree.maxWeight = 1.0
        tree.prodXSec = 1.0*fb
        pvIndex = tree.add_node(pv)
        n1Index = tree.add_node(n1)
        c2Index = tree.add_node(c2)
        tree.add_edge(pvIndex,n1Index)
        tree.add_edge(pvIndex,c2Index)


        # Validate subtree cache for C2 and descendants.
        particleOrderDict = get_particle_order_dict(model)
        decaysDict = get_lightweight_decays(model, particleOrderDict)
        

        pvID = hash(pv.particle)
        c2ID = hash(c2.particle)
        n2ID = hash(n2.particle)
        n1ID = hash(n1.particle)

        pvDecay = decayTupleObj(mom=pvID, daughters=[c2ID,n1ID], br=1.0)
        decaysDict[pvID] = [pvDecay]

        cache = build_subtree_cacheFor(pvID, decaysDict, particleOrderDict,
                                       memo={}, minBR=0.0)

        self.assertIn(c2ID, cache)
        self.assertIn(n2ID, cache)
        self.assertIn(n1ID, cache)

        self.assertEqual(len(cache[pvID]),63)
        self.assertEqual(len(cache[n2ID]),52)
        self.assertEqual(len(cache[n1ID]),1)

        self.assertTrue(all(st.particleIDs[0] == c2ID for st in cache[c2ID]))
        self.assertTrue(all(st.particleIDs[0] == n2ID for st in cache[n2ID]))

        totalBR_c2 = sum(st.decayBRs for st in cache[c2ID])
        totalBR_n2 = sum(st.decayBRs for st in cache[n2ID])
        self.assertTrue(np.isclose(totalBR_c2,1.0,atol=0,rtol=1e-4))
        self.assertTrue(np.isclose(totalBR_n2,1.0,atol=0,rtol=1e-4))

        # Classify top-level decays from subtree cache (root daughters of C2)
        labelDict = {hash(p): str(p) for p in model.SMparticles + model.BSMparticles}
        nN2 = 0
        nC1 = 0
        nN1 = 0
        nOther = 0
        for st in cache[pvID]:
            daughters = [labelDict[pid] for pid in st.particleIDs]
            if 'N2'in daughters:
                nN2 +=1
            elif 'C1+' in daughters or 'C1-' in daughters:
                nC1 +=1 
            elif 'N1'in daughters:
                nN1 +=1 
            else:
                nOther +=1 
        self.assertEqual(nN2,52)
        self.assertEqual(nC1,10)
        self.assertEqual(nN1,1)
        self.assertEqual(nOther,0)

        totalXsec = tree.maxWeight*totalBR_c2*fb
        self.assertAlmostEqual(totalXsec.asNumber(fb),1.0,places=4)

        # Compare with simple count from SLHA file:
        totalDecays = 1    
        for n in tree.nodes:
            if str(n) == 'PV':
                continue
            totalDecays *= getNdecays(model,n.particle.pdg)
        self.assertEqual(totalDecays,len(cache[pvID]))


    def test_decomp_stable(self):

        slhafile = './testFiles/slha/lightEWinos_simple.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-4*GeV)

        pv = ParticleNode(model.getParticlesWith(label='PV')[0])
        c1p = ParticleNode(model.getParticlesWith(pdg=1000024)[0])
        c1m = ParticleNode(model.getParticlesWith(pdg=-1000024)[0])


        tree = TheorySMS()
        tree.maxWeight = 1.0
        tree.prodXSec = 1.0*fb
        pvIndex = tree.add_node(pv)
        c1pIndex = tree.add_node(c1p)
        c1mIndex = tree.add_node(c1m)
        tree.add_edge(pvIndex,c1pIndex)
        tree.add_edge(pvIndex,c1mIndex)

        particleOrderDict = get_particle_order_dict(model)
        decaysDict = get_lightweight_decays(model, particleOrderDict)

        pvID = hash(pv.particle)
        c1pID = hash(c1p.particle)
        c1mID = hash(c1m.particle)

        pvDecay = decayTupleObj(mom=pvID, daughters=[c1pID,c1mID], br=1.0)
        decaysDict[pvID] = [pvDecay]

        cache = build_subtree_cacheFor(pvID, decaysDict, particleOrderDict, memo={}, minBR=0.0)

        treeList = cache[pvID]
        nExpected = 2*2 + 1*2 + 2*1 +1
        self.assertEqual(len(treeList),nExpected)

    def test_decomp_stree(self):

        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)

        pv = ParticleNode(model.getParticlesWith(label='PV')[0])
        n2 = ParticleNode(model.getParticlesWith(pdg=1000023)[0])
        c2 = ParticleNode(model.getParticlesWith(pdg=1000037)[0])

        tree = TheorySMS()
        tree.maxWeight = 3.5
        tree.prodXSec = 3.5*fb
        pvIndex = tree.add_node(pv)
        n2Index = tree.add_node(n2)
        c2Index = tree.add_node(c2)
        tree.add_edge(pvIndex,n2Index)
        tree.add_edge(pvIndex,c2Index)

        particleOrderDict = get_particle_order_dict(model)
        decaysDict = get_lightweight_decays(model, particleOrderDict)
        c2ID = hash(c2.particle)
        n2ID = hash(n2.particle)
        pvID = hash(pv.particle)

        pvDecay = decayTupleObj(mom=pvID, daughters=[n2ID,c2ID], br=1.0)
        decaysDict[pvID] = [pvDecay]

        cache = build_subtree_cacheFor(pvID, decaysDict, particleOrderDict, memo={}, minBR=0.0)

        treeList = cache[pvID]
        nN2tot = 52
        nC2tot = 63
        self.assertEqual(len(treeList),nC2tot*nN2tot)

        totalXsec = 0*fb
        nCut = 0
        cut = 0.07*fb
        xsecCut = 0.0*fb
        for t in treeList:
            w = tree.maxWeight*t.decayBRs*fb
            totalXsec += w
            if w > cut:
                nCut += 1
                xsecCut += w
        self.assertAlmostEqual(totalXsec.asNumber(fb),3.5,places=5)

        # Apply an equivalent cut on combinations of cached subtrees.
        filtered = [t for t in treeList if tree.maxWeight*t.decayBRs > cut.asNumber(fb)]
        self.assertEqual(len(filtered),nCut)

        totalXsec = 0*fb
        for t in filtered:
            totalXsec += tree.maxWeight*t.decayBRs*fb
        self.assertAlmostEqual(totalXsec.asNumber(fb),xsecCut.asNumber(fb),places=5)

if __name__ == "__main__":
    unittest.main()                         