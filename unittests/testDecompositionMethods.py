#!/usr/bin/env python3

"""
.. module:: testDecompositionMethods
   :synopsis: Testing methods used for decomposition

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
import numpy as np
from smodels.decomposition.theorySMS import TheorySMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition.decomposer import getDecayNodes, addOneStepDecays,cascadeDecay
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
        decs = getDecayNodes(n2)
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
        # Convert result from particle nodes to strings:
        decs = sorted([ [str(dec[0])] + sorted([str(d) for d in dec[1]]) + [dec[2]] 
                       for dec in decs[:]])
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


        # Check decays step by step:
        tree1Step = addOneStepDecays(tree)
        self.assertEqual(len(tree1Step),4)
        # Count decays:
        nN2 = 0
        nC1 = 0
        nN1 = 0
        nOther = 0
        for t in tree1Step:
            daughters = [str(d) for d in t.daughters(c2Index)]
            if 'N2'in daughters:
                nN2 +=1
            elif 'C1+' in daughters or 'C1-' in daughters:
                nC1 +=1 
            elif 'N1'in daughters:
                nN1 +=1 
            else:
                nOther +=1 
        self.assertEqual(nN2,1)
        self.assertEqual(nC1,2)
        self.assertEqual(nN1,1)
        self.assertEqual(nOther,0)

        tree2Step = []
        for t in tree1Step:
            tree2Step += addOneStepDecays(t)
        self.assertEqual(len(tree2Step),30)
        # Count decays:
        nN2 = 0
        nC1 = 0
        nN1 = 0
        nOther = 0
        for t in tree2Step:
            daughters = [str(d) for d in t.daughters(c2Index)]
            if 'N2'in daughters:
                nN2 +=1
            elif 'C1+' in daughters or 'C1-' in daughters:
                nC1 +=1 
            elif 'N1'in daughters:
                nN1 +=1 
            else:
                nOther +=1 
        self.assertEqual(nN2,20)
        self.assertEqual(nC1,10)
        self.assertEqual(nN1,0)
        self.assertEqual(nOther,0)


        tree3Step = []
        for t in tree2Step:
            daughters = [str(d) for d in t.daughters(c2Index)]
        #     print(daughters)
            tree3Step += addOneStepDecays(t)
        self.assertEqual(len(tree3Step),40)
        # Count decays:
        nN2 = 0
        nC1 = 0
        nN1 = 0
        nOther = 0
        for t in tree3Step:
            daughters = [str(d) for d in t.daughters(c2Index)]
            if 'N2'in daughters:
                nN2 +=1
            elif 'C1+' in daughters or 'C1-' in daughters:
                nC1 +=1 
            elif 'N1'in daughters:
                nN1 +=1 
            else:
                nOther +=1 
        self.assertEqual(nN2,40)
        self.assertEqual(nC1,0)
        self.assertEqual(nN1,0)
        self.assertEqual(nOther,0)


        # Compute full cascade decay
        treeList = cascadeDecay(tree)
        self.assertEqual(len(treeList),63)
        # Count decays:
        nN2 = 0
        nC1 = 0
        nN1 = 0
        nOther = 0
        for t in treeList:
            daughters = [str(d) for d in t.daughters(c2Index)]
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

        totalXsec = 0*fb
        for t in treeList:
            totalXsec += t.maxWeight*fb
        self.assertAlmostEqual(totalXsec.asNumber(fb),1.0,places=4)

        # Compare with simple count from SLHA file:
        totalDecays = 1    
        for n in tree.nodes:
            if str(n) == 'PV':
                continue
            totalDecays *= getNdecays(model,n.particle.pdg)
        self.assertEqual(totalDecays,len(treeList))


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

        treeList = cascadeDecay(tree)
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

        treeList = cascadeDecay(tree)
        nN2tot = 52
        nC2tot = 63
        self.assertEqual(len(treeList),nC2tot*nN2tot)

        totalXsec = 0*fb
        nCut = 0
        cut = 0.07*fb
        xsecCut = 0.0*fb
        for t in treeList:
            w = t.maxWeight*fb
            totalXsec += w
            if w > cut:
                nCut += 1
                xsecCut += w
        self.assertAlmostEqual(totalXsec.asNumber(fb),3.5,places=5)

        # Now change sigmacut:
        treeList = cascadeDecay(tree,sigmacutFB=cut.asNumber(fb))
        self.assertEqual(len(treeList),nCut)

        totalXsec = 0*fb
        for t in treeList:
            t.weightList = t.computeWeightList()
            totalXsec += t.weightList
        self.assertAlmostEqual(totalXsec.asNumber(fb),xsecCut.asNumber(fb),places=5)

if __name__ == "__main__":
    unittest.main()                         