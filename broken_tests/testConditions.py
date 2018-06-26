#!/usr/bin/env python

"""
.. module:: testConditions
   :synopsis: Tests condition violation evaluation

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from databaseLoader import database

class ConditionTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%sinputFiles/slha/lightEWinos.slha" % (installDirectory() )
        topolist = slhaDecomposer.decompose(filename,doCompress=True, doInvisible=True, minmassgap = 5*GeV)
        analyses = database.getExpResults (txnames=["TChiWZoff"])
        theoryPrediction = theoryPredictionsFor(analyses[0], topolist)[0]
        conditionViolation = theoryPrediction.conditions
        self.assertEqual(conditionViolation['Cgtr([[[mu+,mu-]],[[l,nu]]],[[[e+,e-]],[[l,nu]]])'],0.)
        
if __name__ == "__main__":
    unittest.main()
