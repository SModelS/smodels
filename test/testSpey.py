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
        db = Database ( "mldb/" )
        res = db.getExpResults ( "ATLAS-SUSY-2018-04" )[0]
        slhafile = "testFiles/slha4ml/TStauStau_300_106_300_106.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(slhafile)
        smstoplist = decomposer.decompose(model, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        print ( "smstoplist", smstoplist )
        predictions = theoryPredictionsFor ( res, smstoplist )
        print ( "predictions", predictions )


if __name__ == "__main__":
    unittest.main()
