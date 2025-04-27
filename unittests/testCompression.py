#!/usr/bin/env python3

"""
.. module:: testCompression
   :synopsis: Checks the compression algorithms
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.decomposition import decomposer
from smodels.base.physicsUnits import GeV, fb
import unittest
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from unitTestHelpers import canonNameToVertNumb

class CompressionTest(unittest.TestCase):

    def testInvisiblePositive(self):
        """ test the invisible compression, a positive example """
        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile,ignorePromptQNumbers=['spin','eCharge','colordim'])
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
        model.updateParticles(slhafile,ignorePromptQNumbers=['spin','eCharge','colordim'])
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
        model.updateParticles(slhafile,promptWidth = 1e-12*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        topos = decomposer.decompose(model, .1*fb, True, False, 5.*GeV, 5.*GeV)
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

    def testMassISR(self):

        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile,promptWidth = 1e-12*GeV,
                                ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        # Decompose without allowing for ISR compression
        topos = decomposer.decompose(model, 1000*fb, massCompress=True, invisibleCompress=False, 
                                    minmassgap=10.*GeV,minmassgapISR=0*GeV)
        self.assertEqual(len(topos.getSMSList()),14)
        self.assertFalse(110100 in topos) # Check that PV > A + B is not in topos

        # Decompose without allowing for ISR compression up to 1*GeV
        topos = decomposer.decompose(model, 1000*fb, massCompress=True, invisibleCompress=False, 
                                    minmassgap=10.*GeV,minmassgapISR=1*GeV)
        self.assertEqual(len(topos.getSMSList()),14)
        self.assertFalse(110100 in topos) # Check that PV > A + B is not in topos

        # Decompose without allowing for ISR compression up to 5*GeV
        # (now C2-N1 should be compressed)
        topos = decomposer.decompose(model, 1000*fb, massCompress=True, invisibleCompress=False, 
                                    minmassgap=10.*GeV,minmassgapISR=5*GeV)
        self.assertEqual(len(topos.getSMSList()),15)
        self.assertTrue(110100 in topos) # Check that PV > A + B is in topos
        smsISR = topos[110100][0]
        self.assertEqual(str(smsISR),"(PV > N1/N1~,N1/N1~)")



if __name__ == "__main__":
    unittest.main()
