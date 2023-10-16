#!/usr/bin/env python3

"""
.. module:: testSpey
   :synopsis: Tests for the spey interface
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.databaseObj import Database
from smodels.theory import decomposer
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.speyTools import SpeyComputer
from smodels.theory.theoryPrediction import theoryPredictionsFor

class SpeyTest(unittest.TestCase):

    def xestML(self):
        """ test the ML backend """
        db = Database ( "speydb/" )
        res = db.getExpResults ( "ATLAS-SUSY-2018-04" )[0]
        slhafile = "testFiles/slha4spey/TStauStau_300_106_300_106.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist )
        pr = predictions[0]
        lsm = pr.likelihood(0.)
        lbsm = pr.likelihood(1.)
        self.assertAlmostEqual ( lsm, 1.0735609152601552e-43 )
        # lsm, pyhf: 5.626294389030576e-44
        self.assertAlmostEqual ( lbsm, 3.9401532820495447e-45 )
        # lbsm, pyhf: 6.957205346414768e-46
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

    def xestPyhf(self):
        """ test the pyhf backend """
        db = Database ( "speydb/" )
        res = db.getExpResults ( "ATLAS-SUSY-2018-04" )[0]
        slhafile = "testFiles/slha4ml/TStauStau_300_106_300_106.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist )
        pr = predictions[0]
        lsm = pr.likelihood(0.)
        lbsm = pr.likelihood(1.)
        self.assertAlmostEqual ( lsm, 5.626294389030576e-44 )
        # lsm, pyhf: 5.626294389030576e-44
        self.assertAlmostEqual ( lbsm, 6.957205346414768e-46 )
        # lbsm, pyhf: 6.957205346414768e-46
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

    def testPyhfMultipleJsons(self):
        """ test the pyhf backend """
        db = Database ( "speydb/" )
        res = db.getExpResults ( "ATLAS-SUSY-2018-31" )[0]
        slhafile = "testFiles/slha4spey/T6bbHH_720_170_40_720_170_40.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist )
        pr = predictions[0]
        lsm = pr.likelihood(0.)
        lbsm = pr.likelihood(1.)
        self.assertAlmostEqual ( lsm, 5.626294389030576e-44 )
        # lsm, pyhf: 5.626294389030576e-44
        self.assertAlmostEqual ( lbsm, 6.957205346414768e-46 )
        # lbsm, pyhf: 6.957205346414768e-46
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    unittest.main()
