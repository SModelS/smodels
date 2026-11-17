#!/usr/bin/env python3

"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0, "../")

from smodels.statistics.simplifiedLikelihoods import SLLikelihoodComputer
SLLikelihoodComputer.debug_mode = True
from smodels.matching.theoryPrediction import theoryPredictionsFor, \
     TheoryPredictionsCombiner
from smodels.base.smodelsLogging import setLogLevel
from smodels.decomposition import decomposer
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.experiment.databaseObj import Database
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.statistics.basicStats import observed, apriori, aposteriori
from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer
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
        self.assertAlmostEqual(combiner.muhat(), 1.2058828358516187, 4)
        self.assertAlmostEqual(combiner.nllsm(),12.801668574746357, 4)
        self.assertAlmostEqual(combiner.nll(), 12.23876128425637, 4)
        self.assertAlmostEqual(combiner.nll_min(), 12.223407084934427, 4)
        ulmu = combiner.getUpperLimitOnMu()
        # 16.78997035426023/4.71
        self.assertAlmostEqual(ulmu, 3.7926103695052884, 3)
        ulmu_exp = combiner.getUpperLimitOnMu(evaluationType=apriori)
        self.assertAlmostEqual(ulmu_exp, 2.1431580557869347, 3)

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
        defaultnllSMs, defaultnll_min = {}, {}
        # theta_hat = 0., x = 13.
        # scipy.stats.norm.pdf ( x, 13., 3. ) * scipy.stats.poisson.pmf(14, x)
        # = 0.013575602920029094, so we are actually a little off
        defaultnllSMs["ATLAS-CONF-2013-037:SRtN2"] = 4.2840947051889025

        # theta_hat = 2.87723307, x = 3.7 + theta_hat = 6.57723307
        # scipy.stats.norm.pdf(x, 3.7, 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.007423073728232388
        defaultnllSMs["CMS-SUS-16-050-agg:ar8"] = 4.903162058492391

        # nsig = 1., theta_hat = 0., x = 14.
        # scipy.stats.norm.pdf(x, 14.0, 3.0) * scipy.stats.poisson.pmf(14, x)
        # = 0.014094517457734808
        defaultnll_min["ATLAS-CONF-2013-037:SRtN2"] = 4.261969389997849

        # nsig = 5.3, theta_hat = 0, x = 9.
        # scipy.stats.norm.pdf(x, 9., 2.7948166) * scipy.stats.poisson.pmf(9, x)
        # = 0.01880727876784458
        defaultnll_min["CMS-SUS-16-050-agg:ar8"] = 3.973511315574247
        tpreds = theoryPredictionsFor(database, smstopos,
                                      combinedResults=False, useBestDataset=False)
        for t in tpreds:
            t.computeStatistics()
            dId = t.dataset.dataInfo.dataId
            Id = f"{t.dataset.globalInfo.id}:{dId}"
            # print ( "Id", Id )
            nllsm = t.nllsm()
            # print ( "l(mu_hat)", t.likelihood ( 0.03533022229777052 ) )
            # print ( "theta_hat", t.dataset.theta_hat )
            # print ( "dataset", t.dataset.dataInfo.observedN, t.dataset.dataInfo.expectedBG, t.dataset.dataInfo.bgError )
            nll_min = t.nll_min()
            if False:
                print(f"dataset {Id}: theta_hat {t.dataset.theta_hat[0]:.3f} nllsm {nllsm} nll_min {nll_min}")
            # print ( "[er]", Id, "lsm", lsm, "lmax", lmax )
            self.assertAlmostEqual(nllsm, defaultnllSMs[Id], 5)
            self.assertAlmostEqual(nll_min, defaultnll_min[Id], 5)
        # combination:
        # mu_hat 0.035 lmax 0.00011 ul_mu 0.27
        combiner = TheoryPredictionsCombiner(tpreds)
        combiner.computeStatistics()
        fmh = combiner.statsComputer.get_five_values(evaluationType=observed,
                return_nll = True )
        mu_hat, nll_min = fmh["muhat"], fmh["nll_min"]
        nllsm = combiner.nllsm()
        # print ( "muhat", mu_hat, "lmax", lmax )
        # multiply the previous lsms, 0.013786096355236995 * 0.007423073728232388
        # = 0.00010233520966944002
        self.assertAlmostEqual(nllsm, 9.187256763681294, 4)
        # mu_hat is determined numerically, but its easy to verify graphically,
        # see http://smodels.github.io/test/testTheoryPredCombinations.png
        self.assertAlmostEqual(mu_hat, 0.03533022229777052, 4)
        # lmax must be the product of likelihoods evaluated at mu_hat
        # 0.007672358984439363 * 0.014016921020572387
        # = 0.00010754284992636553
        self.assertAlmostEqual(nll_min, 9.137621066391294, 4)

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
        #mdbpath = 'unittest+unittestextra'
        from databaseLoader import dbpath, dbpathextra
        mdbpath = f"{dbpath}+{dbpathextra}"
        #if "./database" in dbpath:
        #    # seems like we are meant to use local databases
        #    mdbpath =  './database+./database_extra/'
        db = Database( mdbpath )
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
        # IDs that should be selected and the respective evaluationType r-values:
        goodIDs = {
#            "CMS-SUS-16-036": (1.379, "upperLimit"),
            "CMS-SUS-12-024": (0.0004525513371509076, "efficiencyMap"),
            "ATLAS-SUSY-2018-12": (2.294e-3, "efficiencyMap"),
            "ATLAS-SUSY-2019-09": (0.231855657, "combined"),
        }
        # Make sure each ID appears only once:
        selectedIDs = {tp.analysisId() : (tp.getRValue(evaluationType=apriori),tp.dataType())
                        for tp in combiner.theoryPredictions}
        self.assertEqual(sorted(list(selectedIDs.keys())),sorted(list(goodIDs.keys())))
        # Check if the correct predictions were selected:
        for ana in goodIDs:
            diff_rel = abs(goodIDs[ana][0]-selectedIDs[ana][0])/goodIDs[ana][0]
            if abs ( diff_rel ) > 2e-3:
                from smodels.base.smodelsLogging import logger
                logger.error ( f"r-values differ for {ana}: {goodIDs[ana][0]}!={selectedIDs[ana][0]}" )
            self.assertAlmostEqual(diff_rel,0.,2)
            self.assertEqual(goodIDs[ana][1], selectedIDs[ana][1])

        self.assertAlmostEqual(combiner.nllsm() / 60.08627570224895, 1., 2)
        self.assertAlmostEqual(combiner.nll() / 60.34829758771848, 1., 2)
        self.assertAlmostEqual(combiner.nll_min() / 60.08627570224895, 1., 2)
        self.assertAlmostEqual(combiner.getRValue() / 0.26067132943352256, 1., 2)
        self.assertAlmostEqual(combiner.CLs(), 0.5745589222694297, 2 )
        self.assertAlmostEqual(combiner.CLs( evaluationType = apriori ), 0.6370833948782422, 2 )
        self.assertAlmostEqual(combiner.CLs( mu=.5 ), 0.7752652260987847, 2 )


    def testGetCombiner(self):
        setLogLevel ( "fatal" )

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        parfile = "./testParameters_comb.ini"

        combiner = getCombiner(slhafile, parfile, database = database )
        nll_min = combiner.nll_min()
        nllsm = combiner.nllsm()
        nllbsm = combiner.nll(mu=1.0)
        nllbsmE = combiner.nll(mu=1.0, evaluationType=apriori)

        nllbsm = combiner.nll( mu=1.0)
        nllbsmE = combiner.nll( mu=1.0, evaluationType=apriori)

        self.assertAlmostEqual(nllsm, 2.293189355862929, 2)
        self.assertAlmostEqual(nllbsm, 1.8787336722875243, 2)
        self.assertAlmostEqual(nll_min, 1.550852032050606, 2)
        self.assertAlmostEqual(combiner.getRValue(), 0.11910627080455084, 4)

        # Also check if likelihood dict is defined:
        muvals = np.linspace(0.,3.,10)
        llhdDict = combiner.getLlhds(muvals=muvals,normalize=True)

        # Check if keys agree
        keys = sorted(['combined','CMS-SUS-16-050-agg','CMS-SUS-13-012'])
        self.assertEqual(keys,sorted(list(llhdDict.keys())))

        for anaId,llhvals in llhdDict.items():
            self.assertEqual(len(muvals),len(llhvals)) # Just check if size of lists match
            norm = np.sum(llhvals)
            self.assertAlmostEqual(norm,1.,2) # Check normalization


class _DummyTheoryPrediction(object):
    """A theory prediction with a known analytic nll curve."""

    def __init__(self, anaId, nllfunc):
        self._anaId = anaId
        self._nllfunc = nllfunc

    def analysisId(self):
        return self._anaId

    def computeStatistics(self, evaluationType=observed):
        pass

    def nll(self, mu=1.0, evaluationType=observed, asimov=None):
        try:
            mu = mu[0]
        except (TypeError, IndexError):
            pass
        return self._nllfunc(float(mu))


class _DummyCombiner(AnaCombLikelihoodComputer):
    """An analyses combination with analytic nll curves, no setup required."""

    def __init__(self, nllfunc, theoryPredictions=None):
        self.theoryPredictions = theoryPredictions or []
        self._nllfunc = nllfunc

    def nll(self, mu=1.0, evaluationType=observed, asimov=None):
        try:
            mu = mu[0]
        except (TypeError, IndexError):
            pass
        return self._nllfunc(float(mu))


class NLLNormalizationTest(unittest.TestCase):

    def testGetNLLsAndGetLlhdsNormalization(self):
        muvals = np.linspace(0., 3., 10)
        combined = lambda mu: .5 * (mu - 1.)**2
        analysis = lambda mu: .25 * mu**2 + .1
        combiner = _DummyCombiner(combined,
                [_DummyTheoryPrediction("ana1", analysis)])

        raw = combiner.getNLLs(muvals, normalize=False)
        self.assertEqual(sorted(raw.keys()), ["ana1", "combined"])
        np.testing.assert_allclose(raw['combined'], [combined(m) for m in muvals])
        np.testing.assert_allclose(raw['ana1'], [analysis(m) for m in muvals])

        normalized = combiner.getNLLs(muvals, normalize=True)
        for anaID, nll in normalized.items():
            # exp(-nll) sums to one for every entry of the dictionary
            self.assertAlmostEqual(np.sum(np.exp(-nll)), 1., 10)
            # normalization is a constant shift in nll space
            np.testing.assert_allclose(nll - raw[anaID],
                    np.full_like(nll, np.mean(nll - raw[anaID])), atol=1e-10)

        llhds, nlls = combiner.getLlhds(muvals, normalize=True,
                                        return_nlls=True)
        for anaID in normalized:
            np.testing.assert_allclose(llhds[anaID], np.exp(-normalized[anaID]))
            np.testing.assert_allclose(nlls[anaID], normalized[anaID])

        unnormalized = combiner.getLlhds(muvals, normalize=False)
        np.testing.assert_allclose(unnormalized['combined'],
                np.exp(-raw['combined']))

    def testNormalizationWithLargeNLLs(self):
        """A plain log(sum(exp(-nll))) would underflow for these nlls."""
        muvals = np.linspace(0., 3., 10)
        offset = 1000.
        combined = lambda mu: offset + .5 * (mu - 1.)**2
        combiner = _DummyCombiner(combined)

        raw = combiner.getNLLs(muvals, normalize=False)['combined']
        normalized = combiner.getNLLs(muvals, normalize=True)['combined']
        self.assertTrue(np.all(np.isfinite(normalized)))
        self.assertAlmostEqual(np.sum(np.exp(-normalized)), 1., 10)
        np.testing.assert_allclose(normalized - raw,
                np.full_like(raw, np.mean(normalized - raw)))


if __name__ == "__main__":
    unittest.main()
