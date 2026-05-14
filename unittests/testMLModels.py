#!/usr/bin/env python3

"""
.. module:: testMLModels
   :synopsis: Test the nnInterface

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.statistics.nnAdapter import NNAdapter
from smodels.base.smodelsLogging import logger
import warnings

class MLModelsTest(unittest.TestCase):
    def testSimple(self):
        """
        A simple case we test "by hand"
        """
        onnxFile = "testFiles/test.onnx"
        adapter = NNAdapter ( onnxFile, False )
        yields = {}
        regions = [ 'SRhigh_0Jb_cuts', 'SRhigh_0Jc_cuts', 'SRhigh_0Jd_cuts',
            'SRhigh_0Je_cuts', 'SRhigh_0Jf1_cuts', 'SRhigh_0Jf2_cuts',
            'SRhigh_0Jg1_cuts', 'SRhigh_0Jg2_cuts', 'SRhigh_nJa_cuts',
            'SRhigh_nJb_cuts', 'SRhigh_nJc_cuts', 'SRhigh_nJd_cuts',
            'SRhigh_nJe_cuts', 'SRhigh_nJf_cuts', 'SRhigh_nJg_cuts',
            'SRlow_0Jb_cuts', 'SRlow_0Jc_cuts', 'SRlow_0Jd_cuts',
            'SRlow_0Je_cuts', 'SRlow_0Jf1_cuts', 'SRlow_0Jf2_cuts',
            'SRlow_0Jg1_cuts', 'SRlow_0Jg2_cuts', 'SRlow_nJb_cuts',
            'SRlow_nJc_cuts', 'SRlow_nJd_cuts', 'SRlow_nJe_cuts',
            'SRlow_nJf1_cuts', 'SRlow_nJf2_cuts', 'SRlow_nJg1_cuts',
            'SRlow_nJg2_cuts', 'CR_0J_WZ_cuts', 'CR_nJ_WZ_cuts' ]
        for region in regions: # predict for no yields
            yields[ region ] = 0.
        ret = adapter.predict ( yields )
        truths = { 'nll_exp_0': 675.10349848,
                   'nll_exp_1': 675.2284106540345,
                   'nll_obs_0': 688.4482887699999,
                   'nll_obs_1': 691.6445574490972,
                   'nllA_exp_0': 675.10349845,
                   'nllA_exp_1': 675.228398256828,
                   'nllA_obs_0': 674.7278682388554,
                   'nllA_obs_1': 674.427532987276,
                   'nll_obs_max': 682.9210459104119 }
        for name,value in truths.items():
            self.assertAlmostEqual ( ret[name], value )

    def testRunMLModels(self):
        import warnings

        warnings.filterwarnings(
            "ignore",
            category=DeprecationWarning,
            message=r".*RefResolver is deprecated.*",
            module=r"pyhf\.schema\.validator",
        )

        from smodels.experiment.databaseObj import Database
        from smodels.decomposition import decomposer
        from smodels.base import runtime
        from smodels.tools.particlesLoader import load
        from smodels.base.model import Model
        from smodels.base.physicsUnits import GeV
        import os
        from smodels.share.models.SMparticles import SMList
        from smodels.matching.theoryPrediction import theoryPredictionsFor
        from smodels.statistics.basicStats import observed, apriori, aposteriori
        db = Database ( "./mlmodels_db/" ) ## a small database with mlmodels
        db.getExpResults()
        slhafile = os.path.abspath('./testFiles/slha/TChiWZoff_150_125_150_125.slha')
        slhafile = os.path.abspath('./testFiles/mlmodels/ewkinos.slha')
        runtime.modelFile = "smodels.share.models.mssm"
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,
                ignorePromptQNumbers = ['eCharge','colordim','spin'])

        topDict = decomposer.decompose(model, sigmacut=0.001,
                               massCompress=True, invisibleCompress=True,
                               minmassgap=5*GeV)
        allPredictions = theoryPredictionsFor( db, topDict,
                combinedResults=True )
        nlls = { 'ATLAS-SUSY-2019-09': {
            'obs': 109.2548744617272, 'exp': 99.14532271766798 },
                 'ATLAS-SUSY-2018-32': {
            'obs': 93.80734619587368, 'exp': 82.58896814755686 } }
        uls = { 'ATLAS-SUSY-2019-09': {
            'obs': 0.27512461085675616, 'exp': 0.3828992015046798,
            'p1': 0.28588195660979965 },
                 'ATLAS-SUSY-2018-32': {
            'obs': 1.3190277005301172,  'exp': 1.5852880131898859,
            'p1': 1.3302383061447438 } }
        for p in allPredictions:
            if p.nll() != nlls[p.analysisId()]["obs"]:
                print ( f"[testMLModels] for {p.analysisId()}:" )
            self.assertAlmostEqual ( p.nll(),
                    nlls[p.analysisId()]["obs"] )
            self.assertAlmostEqual ( p.nll( evaluationType = apriori),
                    nlls[p.analysisId()]["exp"] )
            self.assertAlmostEqual ( p.getUpperLimitOnMu(),
                    uls[p.analysisId()]["obs"] )
            self.assertAlmostEqual ( \
                    p.getUpperLimitOnMu( evaluationType = aposteriori ),
                    uls[p.analysisId()]["exp"] )
            self.assertAlmostEqual ( p.getUpperLimitOnMu( pmSigma = 1 ),
                    uls[p.analysisId()]["p1"] )

if __name__ == "__main__":
    unittest.main()
