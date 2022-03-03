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


class CombinedTheoryPredsTest(unittest.TestCase):

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [f, f.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)

    def testConstruction(self):
        """ this method should simply test if the fake result and the
            covariance matrix are constructed appropriately """
        database = Database("unittest")
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
        database = Database("unittest")
        dTypes = ["efficiencyMap"]
        anaids = [ "CMS-SUS-13-012", "ATLAS-CONF-2013-037" ]
        dsids = [ "3NJet6_1250HT1500_300MHT450", "SRtN2" ]
        # CMS-SUS-13-012
        # dataId: 3NJet6_1250HT1500_300MHT450
        # observedN: 38
        # expectedBG: 42.8
        # bgError: 9.5
        #
        # ATLAS-CONF-2013-037
        # dataId: SRtN2
        # dataType: efficiencyMap
        # observedN: 14
        # expectedBG: 13.0
        # bgError: 3.0
        slhafile = "testFiles/slha/T1tttt.slha"
        exp_results = database.getExpResults(analysisIDs=anaids,
                                             datasetIDs=dsids, dataTypes=dTypes)
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        smstopos = decomposer.decompose(model)
        tpreds = []
        defaultLSMs, defaultLmax = { }, { }
        # FIXME add the exact derivation
        # mention theta_hat
        # poisson ( 38, x ) * gauss ( x, 42.8, 9.5 )
        defaultLSMs["CMS-SUS-13-012:SRtN2"] = 0.013786096355236995
        defaultLSMs["CMS-SUS-13-012:3NJet6_1250HT1500_300MHT450" ] = 0.0024804685610203808
        defaultLmax["CMS-SUS-13-012:SRtN2"] = 0.014094517457734808
        defaultLmax["CMS-SUS-13-012:3NJet6_1250HT1500_300MHT450" ] = 0.0024804685610203808
        for er in exp_results:
            ts = theoryPredictionsFor(er, smstopos,
                combinedResults=False, useBestDataset=False, marginalize=False)
            for t in ts:
                tpreds.append ( t )
        for t in tpreds:
            t.computeStatistics()
            dId = t.dataset.dataInfo.dataId
            Id = f"{er.globalInfo.id}:{dId}"
            lsm = t.lsm()
            # print ( "dataset", t.dataset.theta_hat )
            lmax = t.lmax()
            # print ( "dataset", t.dataset.theta_hat )
            # print ( "[er]", Id, "lsm", lsm, "lmax", lmax )
            self.assertAlmostEqual ( lsm, defaultLSMs[Id], 5 )
            self.assertAlmostEqual ( lmax, defaultLmax[Id], 5 )

if __name__ == "__main__":
    unittest.main()
