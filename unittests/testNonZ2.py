#!/usr/bin/env python3

"""
.. module:: testNonZ2
   :synopsis: Tests Non-Z2 topologies

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import os
sys.path.insert(0, "../")
import unittest
from smodels.base import runtime
from smodels.tools.particlesLoader import load
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition import decomposer
from smodels.base.physicsUnits import fb, GeV, pb
from smodels.base.smodelsLogging import setLogLevel
from smodels.matching.theoryPrediction import theoryPredictionsFor
from databaseLoader import database
setLogLevel("error")


class RunNonZ2Test(unittest.TestCase):

    def testTRV1(self):

        database.selectExpResults()
        runtime.modelFile = "./testFiles/slha/TRV1_1800_300_300.slha"
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        slhafile = './testFiles/slha/TRV1_1800_300_300.slha'
        model.updateParticles(inputFile=slhafile,ignorePromptQNumbers = ['eCharge','colordim'],promptWidth=1e-5*GeV)

        
        # Set main options for decomposition
        sigmacut = 0.01*fb
        mingap = 5.*GeV

        # Decompose model
        topDict = decomposer.decompose(model, sigmacut,
                                    massCompress=True, invisibleCompress=True,
                                    minmassgap=mingap)
        
        self.assertEqual(len(topDict.getSMSList()),3)

        weights = sorted([0.908*pb, 0.908*pb, 0.39*pb])
        decompW = sorted([sms.weightList[0].value for sms in topDict.getSMSList()])
        for iw,w in enumerate(weights):
            self.assertAlmostEqual(w.asNumber(pb),decompW[iw].asNumber(pb),places=2)

        allPredictions = theoryPredictionsFor(database, topDict, combinedResults=False)

        self.assertEqual(len(allPredictions),1)

        p = allPredictions[0]
        self.assertAlmostEqual(p.xsection.asNumber(pb),0.39,places=2)
        self.assertAlmostEqual(p.getRValue(),1.744,places=2)



if __name__ == "__main__":
    unittest.main()
