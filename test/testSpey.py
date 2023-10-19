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

    def testML(self):
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
        lsm = pr.likelihood(0.,return_nll=True)
        lbsm = pr.likelihood(1.,return_nll=True)
        self.assertAlmostEqual ( lsm, 99.58629305653673, 4 )
        # lsm, pyhf: 5.626294389030576e-44
        self.assertAlmostEqual ( lbsm, 103.97913641499024, 4 )
        # lbsm, pyhf: 6.957205346414768e-46
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

    def testPyhf(self):
        """ test the pyhf backend """
        db = Database ( "speydb/" )
        res = db.getExpResults ( "ATLAS-SUSY-2018-04" )[0]
        slhafile = "testFiles/slha4spey/TStauStau_300_106_300_106.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist, combinedResults=True )
        pr = predictions[0]
        lsm = pr.likelihood(0.,return_nll=True)
        lbsm = pr.likelihood(1.,return_nll=True)
        self.assertAlmostEqual ( lsm, 99.58629305653673, 4 )
        self.assertAlmostEqual ( lbsm, 103.97913641499024, 4 )

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
        lsm = pr.likelihood(0.,return_nll=True)
        lbsm = pr.likelihood(1.,return_nll=True)
        self.assertAlmostEqual ( lsm, 3.76268804500719, 4 )
        # lsm, spey: 0.013157507015802185
        self.assertAlmostEqual ( lbsm, 5.2374904021396596, 4 )
        # lbsm, pyhf: 6.957205346414768e-46
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

    def testCombo(self):
        """ test the analysis combination bit """
        db = Database ( "database" )
        res = db.getExpResults ( "CMS-SUS-13-012" )[0]
        slhafile = "./testFiles/slha/gluino_squarks.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist )
        pr = predictions[0]
        lsm = pr.likelihood(0.)
        lbsm = pr.likelihood(1.)
        # the 5.2 == bgErr is due to the different parametrizations
        self.assertAlmostEqual ( lsm, 0.002530290884739633 * 5.2, 4 )
        # lsm, spey: 0.013157507015802185
        # lsm, SL: 0.002530290884739633
        self.assertAlmostEqual ( lbsm, 0.0038277234526390034 * 5.2, 4 )
        # lbsm, spey: 0.01990416119128329
        # lbsm, SL: 0.0038277234526390034
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()

    def testSLv2(self):
        """ test the pyhf backend """
        db = Database ( "speydb/" )
        res = db.getExpResults ( "CMS-SUS-20-004" )[0]
        slhafile = "testFiles/slha4spey/TChiHH_400_100_400_100.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        predictions = theoryPredictionsFor ( res, smstoplist, combinedResults=True )
        pr = predictions[0]
        lsm = pr.likelihood(0.,return_nll=True)
        lbsm = pr.likelihood(1.,return_nll=True)
        # print ( "predictions", predictions, "lsm", lsm )
        # log(sqrt(det)) = 7.20251460179376
        self.assertAlmostEqual ( lsm, 66.90365456784738 , 4 ) # spey v2
        self.assertAlmostEqual ( lbsm, 63.988519076264666, 4 ) # spey v2
        # delta in spey is 2.915
        #self.assertAlmostEqual ( lsm, 77.80888820981092, 4 ) # SLv2
        # self.assertAlmostEqual ( lbsm, 75.06901967038128, 4 ) # SLv2
        # delta in SLv2 is 2.74
        #self.assertAlmostEqual ( lsm, 67.60981938110237, 4 ) # spey v1
        # self.assertAlmostEqual ( lbsm, 64.24133538273146, 4 ) # spey v1
        # delta in spey v1 is 3.3684



if __name__ == "__main__":
    unittest.main()
