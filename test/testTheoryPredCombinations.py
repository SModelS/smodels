#!/usr/bin/env python3

"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0, "../")

from smodels.tools.simplifiedLikelihoods import LikelihoodComputer
LikelihoodComputer.debug_mode = True
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.theory import decomposer
from smodels.tools.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.theory.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.experiment.databaseObj import Database
from unitTestHelpers import equalObjs, runMain, importModule
import unittest
import os
from databaseLoader import database

class CombinedTheoryPredsTest(unittest.TestCase):

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [f, f.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)

    def testConstruction(self):
        """ this method should simply test if the fake result and the
            covariance matrix are constructed appropriately """
        dTypes = ["efficiencyMap"]
        anaids = ["CMS-SUS-16-050-agg", "CMS-SUS-13-012"]
        dsids = ["ar8", "ar9", "3NJet6_1250HT1500_300MHT450"]
        slhafile = "testFiles/slha/T1tttt.slha"
        exp_results = database.getExpResults(analysisIDs=anaids,
                                             datasetIDs=dsids, dataTypes=dTypes)
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        smstopos = decomposer.decompose(model)
        tpreds = []
        for er in exp_results:
            ts = theoryPredictionsFor(er, smstopos,
                combinedResults=False, useBestDataset=False, marginalize=False)
            for t in ts:
                t.computeStatistics()
                # print("er", str(er), "lsm", t.lsm, "lmax", t.lmax)
                tpreds.append(t)
        combiner = TheoryPredictionsCombiner(tpreds)
        combiner.computeStatistics()
        mu_hat, sigma_mu, lmax = combiner.findMuHat(allowNegativeSignals=True,
                                                    extended_output=True)
        self.assertAlmostEqual(combiner.lsm(), 2.756169857697467e-06, 4)
        self.assertAlmostEqual(combiner.likelihood(), 5.001298746531528e-06, 4)
        self.assertAlmostEqual(combiner.lmax(), 5.131156389020586e-06, 4)
        ulmu = combiner.getUpperLimitOnMu()
        # 16.78997035426023/4.71
        self.assertAlmostEqual(ulmu, 3.45262, 3)

    def testByHandComputed ( self ):
        """ a unit test where in the comments I show the manual computations, step by step, for comparison """
	      # see http://smodels.github.io/test/testTheoryPredCombinations.png
        dTypes = ["efficiencyMap"]
        anaids = [ "CMS-SUS-16-050-agg", "ATLAS-CONF-2013-037" ]
        dsids = [ "SRtN2", "ar8" ]
        # ATLAS-CONF-2013-037
        # dataId: SRtN2
        # dataType: efficiencyMap
        # observedN: 14
        # expectedBG: 13.0
        # bgError: 3.0

        # CMS-SUS-16-050-agg
        # dataId: ar8
        # observedN: 9
        # expectedBG: 3.7
        # bgError: 2.7948166
        slhafile = "testFiles/slha/T1tttt.slha"
        exp_results = database.getExpResults(analysisIDs=anaids,
                                             datasetIDs=dsids, dataTypes=dTypes)
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        smstopos = decomposer.decompose(model)
        tpreds = []
        defaultLSMs, defaultLmax = { }, { }
        # theta_hat = 0., x = 13.
        # scipy.stats.norm.pdf ( x, 13., 3. ) * scipy.stats.poisson.pmf(14, x)
        # = 0.013575602920029094, so we are actually a little off
        defaultLSMs["ATLAS-CONF-2013-037:SRtN2"] = 0.013786096355236995

        # theta_hat = 2.87723307, x = 3.7 + theta_hat = 6.57723307
        # scipy.stats.norm.pdf(x, 3.7, 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.007423073728232388
        defaultLSMs["CMS-SUS-16-050-agg:ar8" ] = 0.007423073728232388

        # nsig = 1., theta_hat = 0., x = 14.
        # scipy.stats.norm.pdf(x, 14.0, 3.0) * scipy.stats.poisson.pmf(14, x)
        # = 0.014094517457734808
        defaultLmax["ATLAS-CONF-2013-037:SRtN2"] = 0.014094517457734808

        # nsig = 5.3, theta_hat = 0, x = 9.
        # scipy.stats.norm.pdf(x, 9., 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.01880727876784458
        defaultLmax["CMS-SUS-16-050-agg:ar8" ] = 0.01880727876784458
        for er in exp_results:
            ts = theoryPredictionsFor(er, smstopos,
                combinedResults=False, useBestDataset=False, marginalize=False)
            for t in ts:
                tpreds.append ( t )
        for t in tpreds:
            t.computeStatistics()
            dId = t.dataset.dataInfo.dataId
            Id = f"{t.dataset.globalInfo.id}:{dId}"
            #print ( "Id", Id )
            lsm = t.lsm()
            # print ( "l(mu_hat)", t.likelihood ( 0.03533022229777052 ) )
            #print ( "theta_hat", t.dataset.theta_hat )
            #print ( "dataset", t.dataset.dataInfo.observedN, t.dataset.dataInfo.expectedBG, t.dataset.dataInfo.bgError )
            lmax = t.lmax()
            if False:
                print ( f"dataset {Id}: theta_hat {t.dataset.theta_hat[0]:.3f} lsm {lsm} lmax {lmax}" )
            # print ( "[er]", Id, "lsm", lsm, "lmax", lmax )
            self.assertAlmostEqual ( lsm, defaultLSMs[Id], 5 )
            self.assertAlmostEqual ( lmax, defaultLmax[Id], 5 )
        # combination:
        # mu_hat 0.035 lmax 0.00011 ul_mu 0.27
        combiner = TheoryPredictionsCombiner(tpreds)
        combiner.computeStatistics()
        mu_hat, sigma_mu, lmax = combiner.findMuHat( allowNegativeSignals=True,
                                                     extended_output=True )
        lsm = combiner.lsm()
        # print ( "muhat", mu_hat, "lmax", lmax )
	      # multiply the previous lsms, 0.013786096355236995 * 0.007423073728232388
	      # = 0.00010233520966944002
        self.assertAlmostEqual ( lsm, 0.00010233520966944002, 4 )
	      # mu_hat is determined numerically, but its easy to verify graphically,
	      # see http://smodels.github.io/test/testTheoryPredCombinations.png
        self.assertAlmostEqual ( mu_hat, 0.03533022229777052, 4 )
	      # lmax must be the product of likelihoods evaluated at mu_hat
	      # 0.007672358984439363 * 0.014016921020572387
	      # = 0.00010754284992636553
        self.assertAlmostEqual ( lmax, 0.00010754284992636553, 4 )

if __name__ == "__main__":
    unittest.main()
