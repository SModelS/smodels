#!/usr/bin/env python3

"""
.. module:: testDecomposer
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition import decomposer
from smodels.base.physicsUnits import fb
from smodels.experiment.databaseObj import Database
from smodels.matching.theoryPrediction import theoryPredictionsFor


class ErasePromptTest(unittest.TestCase):

    def testDecomposition(self):

        # Decompose erase quantum numbers 
        # (prompt decaying particles with the same mass will be considered equal)
        filename = "./testFiles/slha/squarks_degenerate.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,ignorePromptQNumbers=['spin','eCharge','colordim'])
        topList = decomposer.decompose(model, sigmacut=3*fb)
        self.assertEqual(len(topList.getSMSList()),1)

        # Decompose keeping quantum numbers 
        # (prompt decaying particles with the same mass, but distinct quantum numbers
        #  will not be considered equal)
        filename = "./testFiles/slha/squarks_degenerate.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
        topListB = decomposer.decompose(model, sigmacut=3*fb)
        self.assertEqual(len(topListB.getSMSList()),2)

        self.assertEqual(len(topList[11101001101000]),1)
        self.assertEqual(len(topListB[11101001101000]),2)
        smsA = topList[11101001101000][0]
        smsB = topListB[11101001101000][0]


        # Only neutralino and SM particles should have charges:
        chargesKeptA = {'N1': 0, 'q': sorted([-2./3., -1./3., 1./3., 2./3.]),'N1~' : 0}
        # Only squarks with same charge should be combined:
        chargesKeptB = {'sc_L/sc_R/su_L/su_R': 2./3., 'sc_L~/sc_R~/su_L~/su_R~': -2./3., 
                        'N1': 0, 'q': sorted([-2./3., -1./3., 1./3., 2./3.]),
                        'N1~' : 0}
        
        
        eChargesA = {}
        for node in smsA.nodes:
            try:
                eChargesA[str(node)] = node.eCharge
            except AttributeError:
                pass # Should fail for BSM particles
        self.assertEqual(len(chargesKeptA),len(eChargesA))
        # Explicitly check the charges
        for ptc,charges in eChargesA.items():
            if isinstance(charges,list):
                charges = sorted(list(set(charges)))
                for ic,c in enumerate(charges):
                    self.assertAlmostEqual(c,chargesKeptA[ptc][ic],
                                           places=4)
            else:
                self.assertAlmostEqual(charges,chargesKeptB[ptc],
                                       places=4)


        eChargesB = {}
        for node in smsB.nodes:
            try:
                eChargesB[str(node)] = node.eCharge
            except AttributeError:
                pass # Should only fail for PV
        self.assertEqual(len(chargesKeptB),len(eChargesB))
        # Explicitly check the charges
        for ptc,charges in eChargesB.items():
            if isinstance(charges,list):
                charges = sorted(list(set(charges)))
                for ic,c in enumerate(charges):
                    self.assertAlmostEqual(c,chargesKeptB[ptc][ic],
                                           places=4)
            else:
                self.assertAlmostEqual(charges,chargesKeptB[ptc],
                                       places=4)
                

    def testTheoryPreds(self):
    

        db = Database('./dbadd1/')
        # Decompose erase quantum numbers 
        # (prompt decaying particles with the same mass will be considered equal)
        filename = "./testFiles/slha/squarks_degenerate.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,ignorePromptQNumbers=['spin','eCharge','colordim'])
        topListA = decomposer.decompose(model, sigmacut=3*fb)
        allPredictionsA = theoryPredictionsFor(db, topListA, combinedResults=False)
        self.assertEqual(len(allPredictionsA),2)
        

        # Decompose keeping quantum numbers 
        # (prompt decaying particles with the same mass, but distinct quantum numbers
        #  will not be considered equal)
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
        topListB = decomposer.decompose(model, sigmacut=3*fb)
        allPredictionsB = theoryPredictionsFor(db, topListB, combinedResults=False)
        self.assertEqual(len(allPredictionsB),2)
        

        # Test UL
        tpA = allPredictionsA[0]
        tpB = allPredictionsB[0]
        self.assertEqual(len(tpA.smsList),1)
        self.assertEqual(len(tpB.smsList),2)
        # Now make sure both theory predictions and r-values are equal:
        self.assertAlmostEqual(tpA.xsection.asNumber(fb),
                               tpB.xsection.asNumber(fb),places=3)

        self.assertAlmostEqual(tpA.getRValue(),
                               tpB.getRValue(),places=4)
        
        # Test EM
        tpA = allPredictionsA[1]
        tpB = allPredictionsB[1]
        self.assertEqual(len(tpA.smsList),1)
        self.assertEqual(len(tpB.smsList),2)
        # Now make sure both theory predictions and r-values are equal:
        self.assertAlmostEqual(tpA.xsection.asNumber(fb),
                               tpB.xsection.asNumber(fb),places=3)

        self.assertAlmostEqual(tpA.getRValue(),
                               tpB.getRValue(),places=4)        
        
        


if __name__ == "__main__":
    unittest.main()
