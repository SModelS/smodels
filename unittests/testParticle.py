#!/usr/bin/env python3

"""
.. module:: testParticleClass
   :synopsis: Tests the smodels.base.particleClass.Particles and .MultiParticle class
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""
import sys
sys.path.insert(0,"../")
import unittest
import numpy as np
from smodels.base.particle import Particle, MultiParticle
from smodels.base.physicsUnits import GeV
from smodels.base.model import Model
from smodels.experiment.expAuxiliaryFuncs import smsInStr
from smodels.experiment.defaultFinalStates import finalStates
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList



p1 = Particle(isSM=False, label='p1', pdg=None, mass=100.*GeV,
              eCharge=None, colordim=None, spin=None, width=None, branches=None)
p2 = Particle(isSM=False, label='p1', pdg=1000021, mass=50*GeV,
              eCharge=None, colordim=None, spin=None, width=None, branches=None)
p3 = Particle(isSM=False, label='p3', pdg=1, mass=110.*GeV,
              eCharge=None, colordim=None, spin=None, width=None, branches=None)
p4 = Particle(isSM=False, label='p4', pdg=None, mass=110.*GeV,
              eCharge=None, colordim=None, spin=None, width=None, branches=None)
p4a = Particle(isSM=False, label='p4~', pdg=None, mass=110.*GeV,
               eCharge=None, colordim=None, spin=None, width=None, branches=None)
p5 = Particle(isSM=False, label='p5+', pdg=2, mass=110.*GeV,
              eCharge=1., colordim=None, spin=None, width=None, branches=None)
p5m = Particle(isSM=False, label='p5-', pdg=-2, mass=110.*GeV,
               eCharge=-1., colordim=None, spin=None, width=None, branches=None)

p1c = p1.copy()
p1c.pdg = 10


class ParticleTest(unittest.TestCase):

    def testSLHAFileAsString(self):
        """ test if we get the same result if SLHA file is given as string
        """
        slhafile = 'testFiles/slha/lightEWinos.slha'
        firstModel = Model(BSMparticles=BSMList, SMparticles=SMList)
        firstModel.updateParticles(slhafile)

        secondModel = Model(BSMparticles=BSMList, SMparticles=SMList)
        with open(slhafile) as f:
            secondModel.updateParticles(f.read())

        self.assertTrue ( firstModel == secondModel )

        slhafile = 'testFiles/slha/gluino_squarks.slha'
        with open(slhafile) as f:
            secondModel.updateParticles(f.read())

        self.assertTrue ( firstModel != secondModel )

    def testParticleComparison(self):

        slhafile = 'testFiles/slha/lightEWinos.slha'
        bsmModel = Model(BSMparticles=BSMList, SMparticles=SMList)
        bsmModel.updateParticles(slhafile)
        BSMparticles = bsmModel.BSMparticles
        SMparticles = bsmModel.SMparticles
        fStates = finalStates.SMparticles
        allParticles = []
        for p in SMparticles+BSMparticles+fStates:
            if any(p is pB for pB in allParticles):
                continue
            allParticles.append(p)
        allParticles = sorted(allParticles, key = lambda p: p._id)
        allIDs = [p._id for p in allParticles]
        for pid in allIDs:
            self.assertTrue(allIDs.count(pid) == 1)
        compMatrixA = np.zeros((len(allParticles),len(allParticles)))
        compMatrixDefault = np.zeros((len(allParticles),len(allParticles)))
        for i,p1 in enumerate(allParticles):
            for j,p2 in enumerate(allParticles):
                compMatrixDefault[i,j] = p1.cmpProperties(p2)
                compMatrixA[i,j] = p1.__cmp__(p2)

        self.assertTrue(np.array_equal(compMatrixDefault,compMatrixA)) #Check if comparison is correct

        compMatrixB = np.zeros((len(allParticles),len(allParticles)))
        for i,p1 in enumerate(allParticles):
            for j,p2 in enumerate(allParticles): #Compare again (now it should use the stored comparisons)
                if p1 == p2:
                    compMatrixB[i,j] = 0
                elif p1 > p2:
                    compMatrixB[i,j] = 1
                else:
                    compMatrixB[i,j] = -1

        self.assertTrue(np.array_equal(compMatrixB,compMatrixA)) #Check if comparison is correct


    def testParticleList(self):
        l1 = MultiParticle(label='plist', particles=[p1,p2])
        lList = finalStates.getParticlesWith(label='l')[0]

        self.assertEqual( l1.label, 'plist')
        self.assertNotEqual( l1 , lList)
        self.assertTrue(l1 == p1)
        self.assertTrue(l1.pdg == [None,1000021] or l1.pdg == [1000021,None])

    def testChargeConjugation(self):
        p5cc = p5.chargeConjugate()
        p3cc = p3.chargeConjugate()
        self.assertEqual(p5cc.label , 'p5-')
        self.assertNotEqual(p5cc,p5)
        self.assertEqual(p5cc,p5m)
        self.assertEqual(p5cc.eCharge , -1.)
        self.assertEqual(p3cc.pdg, -1)
        self.assertEqual(p4.chargeConjugate(), p4a)

    def testInclusiveParticle(self):
        anything = Particle(label='anything')
        l1 = MultiParticle(label='plist', particles=[p1,p2,p4])
        anyEven = finalStates.getParticle(label='anySM')
        anyOdd = finalStates.getParticle(label='anyBSM')
        lList = finalStates.getParticle(label='l')

        self.assertTrue(isinstance(p1, Particle))
        self.assertTrue(p1 == anything)
        self.assertTrue(anything == p1)
        self.assertTrue(anything == l1)
        self.assertTrue(l1 == anything)
        self.assertTrue(anything == lList)
        self.assertTrue(lList == anything)

        self.assertTrue(anything == anyEven)
        self.assertTrue(anything == anyOdd)
        self.assertFalse(anyOdd == anyEven)
        self.assertFalse(anyEven == anyOdd)

    def testInStr(self):
        instring="[[['t+'],['W-']],[['t+'],['W-']]]+[[['t-'],['W+']],[['t-'],['W+']]]+[[['t+'],['W-']],[['t-'],['W+']]]"
        out= smsInStr(instring)
        out = [x.replace("'","") for x in out]
        self.assertEqual(out, ['[[[t+],[W-]],[[t+],[W-]]]', '[[[t-],[W+]],[[t-],[W+]]]', '[[[t+],[W-]],[[t-],[W+]]]'])


if __name__ == "__main__":
    unittest.main()
