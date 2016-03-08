#!/usr/bin/env python

"""
.. module:: testElementID
   :synopsis: Tests the slha checker

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest, sys
sys.path.insert(1, "../")
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV
from smodels.experiment.databaseObjects import Database
from smodels.theory.theoryPrediction import theoryPredictionsFor

class ElementIdTest(unittest.TestCase):
    def testGoodFile(self):

        listOfIDs = [[28, 29, 30, 31, 24, 25, 26, 27],[23]]
        filename = "%sinputFiles/slha/compressedSpec.slha" % (installDirectory() )
        topoList = slhaDecomposer.decompose(filename,doCompress = True, doInvisible=True, minmassgap = 5*GeV)
        database = Database("database/")
        resultlist = database.getExpResults()
        for res in resultlist:
            theorypredictions = theoryPredictionsFor(res, topoList)
            if not theorypredictions: continue
            for tp in theorypredictions:
                self.assertEquals(tp.IDs,listOfIDs.pop())
            

if __name__ == "__main__":
    unittest.main()
