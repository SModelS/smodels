#!/usr/bin/env python

"""
.. module:: testConditions
   :synopsis: Tests condition violation evaluation

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV
from smodels.experiment import smsHelpers, smsAnalysisFactory
from smodels.theory.theoryPrediction import theoryPredictionFor


class ConditionTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%sinputFiles/slha/lightEWinos.slha" % (installDirectory() )
        topolist = slhaDecomposer.decompose(filename,doCompress=True, doInvisible=True, minmassgap = 5*GeV)
        smsHelpers.base="./database/"
        analyses = smsAnalysisFactory.load(topologies="TChiWZoff")
        theoryPrediction = theoryPredictionFor(analyses[0], topolist)[0]
        conditionViolation = theoryPrediction.conditions
        self.assertEqual(conditionViolation['Cgtr([[[mu+,mu-]],[[l,nu]]],[[[e+,e-]],[[l,nu]]])'],0.)
        
if __name__ == "__main__":
    unittest.main()
