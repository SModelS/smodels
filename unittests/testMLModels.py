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
        truths = { 'nll_exp_1' : 751.1195068359,
                   'nll_obs_0' : 688.4482887700,
                   'nll_obs_1' : 763.6027832031,
                   'nllA_exp_0': 675.1034984500,
                   'nllA_exp_1': 751.1370849609,
                   'nllA_obs_0': 674.7278682391,
                   'nllA_obs_1': 750.5723876953 }
        for name,value in truths.items():
            self.assertAlmostEqual ( ret[name], value )

    def mestRunMLModels(self):
        pass

if __name__ == "__main__":
    unittest.main()
