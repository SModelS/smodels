#!/usr/bin/env python3

"""
.. module:: testElementID
   :synopsis: Tests the slha checker

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV
from databaseLoader import database
from smodels.theory.theoryPrediction import theoryPredictionsFor

class ElementIdTest(unittest.TestCase):
    def testGoodFile(self):

        listOfIDs = {'ATLAS-CONF-2013-037': [31, 32, 33, 34, 27, 28, 29, 30], 
                     'ATLAS-SUSY-2013-05' : [26]}
        filename = "./testFiles/slha/higgsinoStop.slha"
        topoList = slhaDecomposer.decompose(filename,doCompress = True, doInvisible=True, minmassgap = 5*GeV)
        resultlist = database.getExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'])
        for res in resultlist:
            theorypredictions = theoryPredictionsFor(res, topoList)
            if not theorypredictions: continue
            self.assertEqual(len(theorypredictions),1)
            tpIDs = theorypredictions[0].IDs 
            self.assertEqual(sorted(tpIDs),sorted(listOfIDs[res.globalInfo.id]))
            

if __name__ == "__main__":
    unittest.main()
