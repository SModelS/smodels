#!/usr/bin/env python3

"""
.. module:: testElementID
   :synopsis: Tests the slha checker

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory import decomposer
from smodels.tools.physicsUnits import GeV,fb
from databaseLoader import database
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model

class ElementIdTest(unittest.TestCase):
    def testGoodFile(self):

        listOfIDs = {'ATLAS-CONF-2013-037':  [31, 26, 33, 34, 27, 28, 29, 30],
                     'ATLAS-SUSY-2013-05' : [32]}
        filename = "./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,promptWidth = 1e-12*GeV) #Make sure C1 and N2 are treated as prompt

        topoList = decomposer.decompose(model, sigcut= 0.1*fb, doCompress=True, doInvisible=True, minmassgap= 5*GeV)
        resultlist = database.getExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'])
        
        for res in resultlist:
            theorypredictions = theoryPredictionsFor(res, topoList)
            if not theorypredictions: continue
            self.assertEqual(len(theorypredictions),1)
            tpIDs = theorypredictions[0].IDs 
            if 'EXO' in res.globalInfo.id: continue #skip EXO searches
            self.assertEqual(sorted(tpIDs),sorted(listOfIDs[res.globalInfo.id]))            

if __name__ == "__main__":
    unittest.main()
