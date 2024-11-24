#!/usr/bin/env python3

"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0, "../")

from smodels.statistics.simplifiedLikelihoods import LikelihoodComputer
LikelihoodComputer.debug_mode = True
from smodels.matching.theoryPrediction import theoryPredictionsFor, \
     TheoryPredictionsCombiner
from smodels.base.smodelsLogging import setLogLevel
from smodels.decomposition import decomposer
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.experiment.databaseObj import Database
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.matching.modelTester import getCombiner
import numpy as np
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
        database.selectExpResults(analysisIDs=anaids,
                                             datasetIDs=dsids, dataTypes=dTypes)
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        smstopos = decomposer.decompose(model)
        tpreds = []
        tpreds = theoryPredictionsFor(database, smstopos,
                combinedResults=False, useBestDataset=False)
        for t in tpreds:
            t.computeStatistics()
        combiner = TheoryPredictionsCombiner(tpreds)
        combiner.computeStatistics()
        self.assertAlmostEqual(combiner.muhat(), 1.205884272994354, 4)
        self.assertAlmostEqual(combiner.lsm(), 2.756169857697467e-06, 4)
        self.assertAlmostEqual(combiner.likelihood(), 5.001298746531528e-06, 4)
        self.assertAlmostEqual(combiner.lmax(), 5.131156389020586e-06, 4)
        ulmu = combiner.getUpperLimitOnMu()
        # 16.78997035426023/4.71
        self.assertAlmostEqual(ulmu, 3.41744, 3)
        ulmu_exp = combiner.getUpperLimitOnMu(expected=True)
        self.assertAlmostEqual(ulmu_exp, 2.143318, 3)

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
        database.selectExpResults(analysisIDs=anaids,
                                             datasetIDs=dsids, dataTypes=dTypes)
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        smstopos = decomposer.decompose(model)
        tpreds = []
        defaultLSMs, defaultLmax = {}, {}
        # theta_hat = 0., x = 13.
        # scipy.stats.norm.pdf ( x, 13., 3. ) * scipy.stats.poisson.pmf(14, x)
        # = 0.013575602920029094, so we are actually a little off
        defaultLSMs["ATLAS-CONF-2013-037:SRtN2"] = 0.013786096355236995

        # theta_hat = 2.87723307, x = 3.7 + theta_hat = 6.57723307
        # scipy.stats.norm.pdf(x, 3.7, 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.007423073728232388
        defaultLSMs["CMS-SUS-16-050-agg:ar8"] = 0.007423073728232388

        # nsig = 1., theta_hat = 0., x = 14.
        # scipy.stats.norm.pdf(x, 14.0, 3.0) * scipy.stats.poisson.pmf(14, x)
        # = 0.014094517457734808
        defaultLmax["ATLAS-CONF-2013-037:SRtN2"] = 0.014094517457734808

        # nsig = 5.3, theta_hat = 0, x = 9.
        # scipy.stats.norm.pdf(x, 9., 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.01880727876784458
        defaultLmax["CMS-SUS-16-050-agg:ar8"] = 0.01880727876784458
        tpreds = theoryPredictionsFor(database, smstopos, 
                                      combinedResults=False, useBestDataset=False)
        for t in tpreds:
            t.computeStatistics()
            dId = t.dataset.dataInfo.dataId
            Id = f"{t.dataset.globalInfo.id}:{dId}"
            # print ( "Id", Id )
            lsm = t.lsm()
            # print ( "l(mu_hat)", t.likelihood ( 0.03533022229777052 ) )
            # print ( "theta_hat", t.dataset.theta_hat )
            # print ( "dataset", t.dataset.dataInfo.observedN, t.dataset.dataInfo.expectedBG, t.dataset.dataInfo.bgError )
            lmax = t.lmax()
            if False:
                print(f"dataset {Id}: theta_hat {t.dataset.theta_hat[0]:.3f} lsm {lsm} lmax {lmax}")
            # print ( "[er]", Id, "lsm", lsm, "lmax", lmax )
            self.assertAlmostEqual(lsm, defaultLSMs[Id], 5)
            self.assertAlmostEqual(lmax, defaultLmax[Id], 5)
        # combination:
        # mu_hat 0.035 lmax 0.00011 ul_mu 0.27
        combiner = TheoryPredictionsCombiner(tpreds)
        combiner.computeStatistics()        
        fmh = combiner.statsComputer.get_five_values(expected=False)
        mu_hat, lmax = fmh["muhat"], fmh["lmax"]
        lsm = combiner.lsm()
        # print ( "muhat", mu_hat, "lmax", lmax )
        # multiply the previous lsms, 0.013786096355236995 * 0.007423073728232388
        # = 0.00010233520966944002
        self.assertAlmostEqual(lsm, 0.00010233520966944002, 4)
        # mu_hat is determined numerically, but its easy to verify graphically,
        # see http://smodels.github.io/test/testTheoryPredCombinations.png
        self.assertAlmostEqual(mu_hat, 0.03533022229777052, 4)
        # lmax must be the product of likelihoods evaluated at mu_hat
        # 0.007672358984439363 * 0.014016921020572387
        # = 0.00010754284992636553
        self.assertAlmostEqual(lmax, 0.00010754284992636553, 4)

    def testFilter(self):
        import warnings
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        warnings.filterwarnings("ignore", category=UserWarning)        
        from smodels.base import runtime
        runtime._experimental["truncatedgaussians"] = True

        anaids = [
          # 'CMS-SUS-16-036',
          'CMS-SUS-12-024',
          'CMS-SUS-12-028',
          'ATLAS-SUSY-2018-12',
          'ATLAS-SUSY-2016-15','ATLAS-SUSY-2019-09']
        db = Database('unittest+unittestextra')
        slhafile = "testFiles/slha/gluino_squarks.slha"
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile)
        sigmacut = 0.005*fb
        mingap = 5.*GeV
        smstopos = decomposer.decompose(model,
                                        sigmacut, massCompress=True, invisibleCompress=True,
                                        minmassgap=mingap)
        tpreds =  theoryPredictionsFor(db, smstopos, 
                                       combinedResults=True, useBestDataset=False)
        combiner = TheoryPredictionsCombiner.selectResultsFrom(tpreds, anaids)
        # IDs that should be selected and the respective expected r-values:
        goodIDs = {
#            "CMS-SUS-16-036": (1.379, "upperLimit"),
            "CMS-SUS-12-024": (4.52551e-4, "efficiencyMap"),
            "ATLAS-SUSY-2018-12": (2.294e-3, "efficiencyMap"),
            "ATLAS-SUSY-2019-09": (0.231855657, "combined"),
        }
        # Make sure each ID appears only once:
        selectedIDs = {tp.analysisId() : (tp.getRValue(expected=True),tp.dataType())
                        for tp in combiner.theoryPredictions}
        self.assertEqual(sorted(list(selectedIDs.keys())),sorted(list(goodIDs.keys())))
        # Check if the correct predictions were selected:
        for ana in goodIDs:
            diff_rel = abs(goodIDs[ana][0]-selectedIDs[ana][0])/goodIDs[ana][0]
            if abs ( diff_rel ) > 1e-3:
                from smodels.base.smodelsLogging import logger
                logger.error ( f"r-values differ for {ana}: {goodIDs[ana][0]}!={selectedIDs[ana][0]}" )
            self.assertAlmostEqual(diff_rel,0.,2)
            self.assertEqual(goodIDs[ana][1], selectedIDs[ana][1])

        self.assertAlmostEqual(combiner.lsm() / 8.032708820262497e-27, 1., 2)
        self.assertAlmostEqual(combiner.likelihood() / 6.1811227091720504e-27, 1., 2)
        self.assertAlmostEqual(combiner.lmax() / 8.032708820262498e-27, 1., 2)
        self.assertAlmostEqual(combiner.getRValue() / 0.2771301968, 1., 2)        
        self.assertAlmostEqual(combiner.CLs(), 0.46721329662185557, 2 )
        self.assertAlmostEqual(combiner.CLs( expected = True ), 0.5295734024482494, 2 )
        self.assertAlmostEqual(combiner.CLs( mu=.5 ), 0.6474406370750514, 2 )


    def testGetCombiner(self):
        setLogLevel ( "fatal" )

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        parfile = "./testParameters_comb.ini"

        combiner = getCombiner(slhafile, parfile)
        lmax = combiner.lmax()
        lsm = combiner.lsm()
        lbsm = combiner.likelihood(mu=1.0)


        self.assertAlmostEqual(lsm, 0.1009, 2)
        self.assertAlmostEqual(lbsm,0.15356, 2)
        self.assertAlmostEqual(lmax, 0.21210, 2)
        self.assertAlmostEqual(combiner.getRValue(), 0.1342, 2)

        # Also check if likelihood dict is defined:
        muvals = np.linspace(0.,3.,10)
        dmu = muvals[1]-muvals[0]
        llhdDict = combiner.getLlhds(muvals,normalize=True)

        # Check if keys agree
        keys = sorted(['combined','CMS-SUS-16-050-agg','CMS-SUS-13-012'])
        self.assertEqual(keys,sorted(list(llhdDict.keys())))

        for anaId,llhvals in llhdDict.items():
            self.assertEqual(len(muvals),len(llhvals)) # Just check if size of lists match
            norm = np.sum(llhvals*dmu)
            self.assertAlmostEqual(norm,1.,2) # Check normalization

if __name__ == "__main__":
    unittest.main()
