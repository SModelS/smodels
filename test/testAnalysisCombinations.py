#!/usr/bin/env python3
 
"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses
 
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
 
"""
 
import sys,os
import importlib
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import equalObjs, runMain, importModule
from smodels.tools.smodelsLogging import logger, setLogLevel
 
class CombinedAnalysisTest(unittest.TestCase):

    def testConstruction(self):
        from smodels.experiment.databaseObj import Database
        database = Database( "unittest" )
        anas_and_sr = { "CMS-PAS-EXO-16-036": [ "c000", "c100" ], 
                        "CMS-SUS-13-012": [ "3NJet6_500HT800_200MHT300" ] }
        anaids = anas_and_sr.keys()
        dsIDs = []
        for anaid,srs in anas_and_sr.items():
            dsIDs += srs
        dsIDs = set ( dsIDs )

        exp_results = database.getExpResults( analysisIDs = anaids,
                                              dataTypes = [ "efficiencyMap" ],
                                              datasetIDs = dsIDs )
        from smodels.tools.simplifiedLikelihoods import SignalRegionsCombiner
        combiner = SignalRegionsCombiner()
        combiner.fromExpResults ( exp_results, anas_and_sr, 
                                  corrs = { "c000": { "c100": .1 } } )
        cov01 = combiner.fakeResult.globalInfo.covariance[0][1]
        self.assertAlmostEqual ( cov01, 6.731e-3 )
        # print ( exp_results )
        #print ( combiner.fakeResult )
        #print ( combiner.covariance )
 
if __name__ == "__main__":
    setLogLevel ( "debug" )
    unittest.main()
