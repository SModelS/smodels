#!/usr/bin/env python

"""
.. module:: testTheoryPrediction
   :synopsis: Testing the theory prediction.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.physicsUnits import fb, GeV, pb
from databaseLoader import database
import inspect
import os

class IntegrationTest(unittest.TestCase):
    def configureLogger(self):
        import logging.config
        fc= inspect.getabsfile(self.configureLogger).replace ( 
                "testTheoryPrediction.py", "integration.conf" )
        logging.config.fileConfig( fname=fc, disable_existing_loggers=False )

    def predictions(self):
        return { 'ATLAS-SUSY-2013-02:T1': 572.168935 * fb }

    def checkAnalysis(self,expresult,smstoplist):
        from smodels.theory.theoryPrediction import theoryPredictionsFor
        theorypredictions = theoryPredictionsFor(expresult, smstoplist)
        defpreds=self.predictions()
        if not theorypredictions:
            print "no theory predictions for",expresult,"??"
            import sys
            sys.exit(-1)
        for pred in theorypredictions:
            m0=str ( int ( pred.mass[0][0]/GeV )  )
            w=pred.value.getDictionary()[(None,None)]
            predval=w['8 TeV (NLL)']
            defpredval=defpreds[expresult.getValuesFor('id')[0]+\
                ":"+pred.txnames[0].getInfo("txName") ]
            self.assertAlmostEqual ( predval / fb, defpredval / fb )

    def testIntegration(self):
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import fb, GeV
        from smodels.theory import slhaDecomposer
        from smodels.experiment.databaseObj import Database
        slhafile = '../inputFiles/slha/simplyGluino.slha'
        self.configureLogger()
        smstoplist = slhaDecomposer.decompose(slhafile, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        listofanalyses = database.getExpResults( 
                analysisIDs= [ "ATLAS-SUSY-2013-02" ], txnames = [ "T1" ] )
        if type(listofanalyses) != list:
            listofanalyses= [ listofanalyses] 
        for analysis in listofanalyses:
            self.checkAnalysis(analysis,smstoplist)

if __name__ == "__main__":
    unittest.main()
