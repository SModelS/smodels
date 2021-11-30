#!/usr/bin/env python3
 
"""
.. module:: testAnalysisCombinations
   :synopsis: Tests the combination of SRs between analyses
 
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
 
"""
 
import unittest
import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.tools.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.theory import decomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor

class CombinedTheoryPredsTest(unittest.TestCase):

    def testConstruction(self):
        """ this method should simply test if the fake result and the 
            covariance matrix are constructed appropriately """
        database = Database( "unittest" )
        dTypes = [ "efficiencyMap" ]
        anaids = [ "CMS-SUS-16-050-agg", "CMS-SUS-13-012" ]
        dsids = [ "ar8", "ar9", "3NJet6_1250HT1500_300MHT450" ]
        slhafile = "testFiles/slha/T1tttt.slha"
        exp_results = database.getExpResults( analysisIDs=anaids, 
                datasetIDs = dsids, dataTypes = dTypes )
        model = Model(BSMparticles=BSMList, SMparticles=SMList )
        model.updateParticles(inputFile=slhafile )
        smstopos = decomposer.decompose ( model )
        tpreds = []
        for er in exp_results:
            ts = theoryPredictionsFor ( er, smstopos, 
                combinedResults=False, useBestDataset=False, marginalize = False )
            for t in ts: 
                t.computeStatistics()
                # print ( "er", str(er), "lsm", t.lsm, "lmax", t.lmax )
                tpreds.append ( t )
        combiner = TheoryPredictionsCombiner ( tpreds )
        combiner.computeStatistics()
        ulmu = combiner.getUpperLimitOnMu()
        self.assertAlmostEqual ( combiner.lsm, 2.756169857697467e-06, 4 )
        self.assertAlmostEqual ( combiner.likelihood, 5.001298746531528e-06, 4 )
        self.assertAlmostEqual ( combiner.lmax, 5.131156389020586e-06, 4 )
        self.assertAlmostEqual ( ulmu, 3.606980755048784, 3 )


if __name__ == "__main__":
    unittest.main()
