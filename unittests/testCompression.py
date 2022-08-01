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
from unitTestHelpers import canonNameToVertNumb

class CompressionTest(unittest.TestCase):

    def testInvisiblePositive(self):
        """ test the invisible compression, a positive example """
        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        topos = decomposer.decompose ( model, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            vertnumb = canonNameToVertNumb(topos,topo)
            if vertnumb !="[][]":
                continue
            for sms in topos[topo]:
                evenParticles = sms.treeToBrackets()[0]
                evenParticles = str(evenParticles).replace("'","").replace(' ', '')
                if evenParticles != "[[],[]]":
                    continue
                tested = True
                trueMothers = [mother for mother in sms.ancestors if mother is not sms]
                if not trueMothers:
                    continue
                evenParticles = trueMothers[0].treeToBrackets()[0]
                evenParticles = str(evenParticles).replace("'","").replace(' ', '')
                self.assertEqual(evenParticles,"[[],[[nu,nu]]]")
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
            vertnumb = canonNameToVertNumb(topos,topo)
            if vertnumb !="[1,1][1,1]":
                continue
            for sms in topos[topo]:
                evenParticles = sms.treeToBrackets()[0]
                evenParticles = str(evenParticles).replace("'","").replace(' ', '')
                if evenParticles !="[[[t+],[t-]],[[q],[W+]]]":
                    continue
                tested = True
                trueMothers = [mother for mother in sms.ancestors if mother is not sms]
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
            vertnumb = canonNameToVertNumb(topos,topo)
            if vertnumb !="[1][1]":
                continue
            for sms in topos[topo]:
                evenParticles = sms.treeToBrackets()[0]
                evenParticles = str(evenParticles).replace("'","").replace(' ', '')
                if evenParticles !="[[[b]],[[b]]]":
                    continue
                mother = sms.ancestors[0]
                masses = [node.mass for node in mother.nodes if node.isSM is False]
                tested = True
                dm=abs(masses[2]-masses[4])/GeV
                #If intermediate BSM states are compared there are two elements ([[[b],[c,q]],[[b],[q,q]]])
                # which do not get combined because their branches differ by the charges of the intermediate states
                self.assertEqual(len(sms.ancestors),24)
                self.assertTrue(dm < 5.0)
        self.assertTrue(tested)

if __name__ == "__main__":
    unittest.main()
