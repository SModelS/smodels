#!/usr/bin/env python3

"""
.. module:: testCompression
   :synopsis: Checks the compression algorithms
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.theory import decomposer
from smodels.tools.physicsUnits import GeV, fb
import unittest
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model


class CompressionTest(unittest.TestCase):

    def testInvisiblePositive(self):
        """ test the invisible compression, a positive example """
        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        topos = decomposer.decompose ( model, .1*fb, False, True, 5.*GeV )
        tested = False
        for topo in topos:
            if str(topo)!="[][]":
                continue
            for element in topo.elementList:
                if str(element)!="[[],[]]":
                    continue
                tested = True
                trueMothers = [mother for mother in element.motherElements if not mother is element]
                if not trueMothers: continue
                self.assertEqual(str(trueMothers[0]),"[[],[[nu,nu]]]")
                self.assertEqual(len(trueMothers), 1)
        self.assertTrue(tested)

    def testInvisibleNegative(self):
        """ test the invisible compression, a negative example """
        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        topos = decomposer.decompose(model, .1*fb, False, True, 5.*GeV)
        tested = False
        for topo in topos:
            if str(topo)!="[1,1][1,1]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[t+],[t-]],[[q],[W+]]]":
                    continue
                tested = True
                trueMothers = [mother for mother in element.motherElements if not mother is element]
                self.assertEqual(len(trueMothers),0)
        self.assertTrue(tested)

    def testMass(self):
        """ test the mass compression, a positive example """
        tested = False
        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile,promptWidth = 1e-12*GeV)
        topos = decomposer.decompose(model, .1*fb, True, False, 5.*GeV)
        for topo in topos:
            if str(topo)!="[1][1]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[b]],[[b]]]":
                    continue
                masses=element.motherElements[0].mass
                tested = True
                dm=abs(masses[0][1]-masses[0][2])/GeV
                #If intermediate BSM states are compared there are two elements ([[[b],[c,q]],[[b],[q,q]]])
                # which do not get combined because their branches differ by the charges of the intermediate states
                self.assertEqual(len(element.motherElements),24)
                self.assertTrue(dm < 5.0)
        self.assertTrue(tested)

if __name__ == "__main__":
    unittest.main()
