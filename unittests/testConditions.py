#!/usr/bin/env python3

"""
.. module:: testConditions
   :synopsis: Tests condition violation evaluation

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.decomposition import decomposer
from smodels.base.model import Model
from smodels.base.physicsUnits import GeV, fb
from smodels.matching.theoryPrediction import theoryPredictionsFor
from databaseLoader import database

class ConditionTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
                
        topolist = decomposer.decompose(model, sigmacut= 0.1*fb, massCompress=True, invisibleCompress=True, minmassgap = 5*GeV)
        database.selectExpResults(txnames=["TChiWZoff"],analysisIDs='ATLAS-SUSY-2013-12')
        theoryPrediction = theoryPredictionsFor(database, topolist)[0]
        conditionViolation = theoryPrediction.conditions
        self.assertAlmostEqual(conditionViolation,[0.],3)
        
if __name__ == "__main__":
    unittest.main()
