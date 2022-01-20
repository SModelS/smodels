#!/usr/bin/env python3

"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0, "../")

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
        ulmu = combiner.getUpperLimitOnMu()
        self.assertAlmostEqual(combiner.lsm(), 2.756169857697467e-06, 4)
        self.assertAlmostEqual(combiner.likelihood(), 5.001298746531528e-06, 4)
        self.assertAlmostEqual(combiner.lmax(), 5.131156389020586e-06, 4)
        self.assertAlmostEqual(ulmu, 16.78997035426023/4.71, 3)

    def testAPI(self):
        
        filename = "testFiles/slha/T1tttt.slha"
        outputfile = runMain(filename, inifile='testParameters_comb.ini', suppressStdout=True)
        smodelsOutput = importModule(outputfile)
        from T1tttt_comb_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'Element', 'database version',
                        'Total missed xsec',
                        'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                        'Total xsec for missing topologies (fb)', 'Total xsec for missing topologies with displaced decays (fb)',
                        'Total xsec for missing topologies with prompt decays (fb)',
                        'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['CombinedRes'] = sorted(smodelsOutputDefault['CombinedRes'],
                                                     key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedDiff=0.02,
                           ignore=ignoreFields, fname=outputfile)
        self.assertTrue(equals)

        for i in ['./output.py', './output.pyc']:
            if os.path.exists(i):
                os.remove(i)
        self.removeOutputs(outputfile)


if __name__ == "__main__":
    unittest.main()
