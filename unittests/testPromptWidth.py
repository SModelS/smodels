#!/usr/bin/env python3

"""
.. module:: testPromptWidth
   :synopsis: Tests width dependence on prompt decays

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
setLogLevel("fatal")

class RunPromptWidth(unittest.TestCase):

    def testWidthInfo(self):

        database.selectExpResults()
        runtime.modelFile = 'mssm'
        filename = "./testFiles/slha/higgsinoStop.slha"
        BSMList = load()
        model = Model(BSMList,SMList)
        model.updateParticles(filename,promptWidth=1e-12*GeV,
                            ignorePromptQNumbers=['spin','eCharge','colordim']) #Force charginos/neutralinos to be considered as prompt
        
        p = model.getParticle(pdg = 1000023)
        self.assertTrue(p.isPrompt())
        self.assertAlmostEqual(p.totalwidth.asNumber(GeV)/1e-9,7.9951,places=2)


    def testWidthInterpolation(self):

        database.selectExpResults()
        runtime.modelFile = "./testFiles/slha/TRV1_1800_300_300.slha"
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        slhafile = './testFiles/slha/TRV1_1800_300_300.slha'
        model.updateParticles(inputFile=slhafile,ignorePromptQNumbers = ['eCharge','colordim'],promptWidth=1e-5*GeV)

        # Set Zprime width
        Zprime = model.getParticle(pdg=55)
        Zprime.totalwidth = 1.780000e+01*GeV

        
        # Set main options for decomposition
        sigmacut = 0.01*fb
        mingap = 5.*GeV

        # Decompose model
        topDict = decomposer.decompose(model, sigmacut,
                                    massCompress=True, invisibleCompress=True,
                                    minmassgap=mingap)
        
        allPredictions = theoryPredictionsFor(database, topDict, combinedResults=False)

        self.assertEqual(len(allPredictions),1)
        p = allPredictions[0]
        self.assertAlmostEqual(p.xsection.asNumber(pb),0.39,places=2)
        self.assertAlmostEqual(p.getRValue(),1.744,places=2)

        
        # Change Zprime width
        Zprime.totalwidth = 5.780000e+01*GeV
        # Decompose model again
        topDict = decomposer.decompose(model, sigmacut,
                                    massCompress=True, invisibleCompress=True,
                                    minmassgap=mingap)
        
        allPredictions = theoryPredictionsFor(database, topDict, combinedResults=False)

        self.assertEqual(len(allPredictions),1)
        p = allPredictions[0]
        self.assertAlmostEqual(p.xsection.asNumber(pb),0.39,places=2)
        self.assertAlmostEqual(p.getRValue(),1.063,places=2)



if __name__ == "__main__":
    unittest.main()
