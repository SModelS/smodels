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
from smodels.tools.signalRegionsCombiner import SignalRegionsCombiner
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
 
class CombinedAnalysisTest(unittest.TestCase):

    def testConstruction(self):
        from smodels.experiment.databaseObj import Database
        database = Database( "unittest" )
        exp_results = database.getExpResults( dataTypes = [ "efficiencyMap" ] )
        combiner = SignalRegionsCombiner()
        labels = [ "CMS-PAS-EXO-16-036:c000", "c100", "CMS-SUS-13-012:3NJet6_500HT800_200MHT300" ]
        datasets = combiner.selectDatasetsFrom ( exp_results, labels )
        #print ( "labels", datasets.keys() )
        combiner.fromDatasets (datasets, corrs = { "c000": { "CMS-PAS-EXO-16-036:c100": .1 } } )
        cov01 = combiner.fakeResult.globalInfo.covariance[1][2]
        # print ( exp_results )
        #print ( combiner.fakeResult )
        #print ( combiner.covariance )
        self.assertAlmostEqual ( cov01, 6.731e-3 )


    def runException ( self ):
        from smodels.experiment.databaseObj import Database
        database = Database( "unittest" )
        exp_results = database.getExpResults( dataTypes = [ "efficiencyMap" ] )
        combiner = SignalRegionsCombiner()
        datasets = combiner.selectDatasetsFrom ( exp_results )
        combiner.fromDatasets (datasets, corrs = { "c000": { "c100": .1 } } )
        cov01 = combiner.fakeResult.globalInfo.covariance[0][1]

    def testException ( self ):
        self.assertRaises ( SModelSError, self.runException )

        # print ( exp_results )
        # print ( combiner.fakeResult )
        # print ( combiner.covariance )
 
if __name__ == "__main__":
    setLogLevel ( "debug" )
    unittest.main()
