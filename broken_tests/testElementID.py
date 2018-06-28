#!/usr/bin/env python

"""
.. module:: testElementID
   :synopsis: Tests the slha checker

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory
from smodels.theory import decomposer
from smodels.tools.physicsUnits import GeV
from databaseLoader import database
from smodels.theory.theoryPrediction import theoryPredictionsFor

class ElementIdTest(unittest.TestCase):
    def testGoodFile(self):

        listOfIDs = {'ATLAS-CONF-2013-037': [28, 29, 30, 31, 24, 25, 26, 27], 
                     'ATLAS-SUSY-2013-05' : [23]}
        filename = "%sinputFiles/slha/higgsinoStop.slha" % (installDirectory() )
        topoList = decomposer.decompose(filename,doCompress = True, doInvisible=True, minmassgap = 5*GeV)
        resultlist = database.getExpResults()
        for res in resultlist:
            theorypredictions = theoryPredictionsFor(res, topoList)
            if not theorypredictions: continue
            self.assertEquals(len(theorypredictions),1)
            tpIDs = theorypredictions[0].IDs
            self.assertEquals(sorted(tpIDs),sorted(listOfIDs[res.globalInfo.id]))
            

if __name__ == "__main__":
    unittest.main()
