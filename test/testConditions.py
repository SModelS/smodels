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
from smodels.theory import decomposer
from smodels.theory.model import Model
from smodels.tools.physicsUnits import GeV, fb
from smodels.theory.theoryPrediction import theoryPredictionsFor
from databaseLoader import database

class ConditionTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "./testFiles/slha/lightEWinos.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)
                
        topolist = decomposer.decompose(model, sigmacut= 0.1*fb, doCompress=True, doInvisible=True, minmassgap = 5*GeV)
        analyses = database.getExpResults(txnames=["TChiWZoff"],analysisIDs='ATLAS-SUSY-2013-12')
        theoryPrediction = theoryPredictionsFor(analyses[0], topolist)[0]
        conditionViolation = theoryPrediction.conditions
        self.assertEqual(conditionViolation['Cgtr([[[mu+,mu-]],[[l,nu]]],[[[e+,e-]],[[l,nu]]])'],0.)
        
if __name__ == "__main__":
    unittest.main()
