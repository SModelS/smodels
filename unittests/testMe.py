#!/usr/bin/env python3

"""
.. module:: testMLModels
   :synopsis: Test the nnInterface

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest

class MeTest(unittest.TestCase):

    def restMe(self):
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
        db = Database ( "./cls_db/" ) ## a small database with mlmodels
        db.getExpResults()
        slhafile = os.path.abspath('../TChiWZ_600_20_600_20.slha' )
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
        for p in allPredictions:
            print ( f"[testMe] for {p.analysisId()}:" )
            print ( "CLs", p.CLs() )
        import sys; import IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    unittest.main()
