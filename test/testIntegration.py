#!/usr/bin/env python

"""
.. module:: testIntegration
   :synopsis: Integration test, tests a simple but complete use case.
              Uses the database in smodels/test/database.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
from smodels.tools.physicsUnits import fb, GeV
import inspect
import os

class IntegrationTest(unittest.TestCase):
    def configureLogger(self):
        import logging.config
        fc= inspect.getabsfile(self.configureLogger).replace ( "testIntegration.py", "integration.conf" )
        logging.config.fileConfig( fname=fc, disable_existing_loggers=False )

    def predictions(self):
        return { 'SUS13006:TChiWZ': 2.42154060952*fb,
                 'SUS12028:T2': None,
                 'SUS12028:T1': 8.43114195825*fb,
                 'SUS12022:TChiChipmSlepL371': 8.43114195825*fb }

    def checkAnalysis(self,analysis,smstoplist):
        from smodels.theory.theoryPrediction import theoryPredictionFor
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        defpreds=self.predictions()
        print "ana",analysis,theorypredictions
        if not theorypredictions:
            return
        ### print(">>>> Ana %s" % analysis.label)
        for pred in theorypredictions:
            # print ( "Pred ana",pred.analysis.label)
            m0=str ( int ( pred.mass[0][0]/GeV )  )
            # print ( "Pred mass0",m0 )
            #print ( "Pred mass0",int(pred.mass[0]/GeV))
            print ( "Pred value",pred.value.getDictionary() )
            predval=pred.value.getDictionary()[(None,None)]['8 TeV (NLL)']
            defpredval=defpreds[pred.analysis.label+m0] 
            # print ( "Pred default value",defpredval )
            self.assertAlmostEqual ( predval / fb, defpredval / fb )

            # print ( "Pred condition",pred.conditions)
            ## print ( "Pred mass",pred.mass)

    def testIntegration(self):
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import fb, GeV
        from smodels.theory import slhaDecomposer
        from smodels.experiment import smsAnalysisFactory, smsHelpers
        smsHelpers.base = installDirectory() + 'test/database/'
        smsHelpers.runs = [ "2012" ]
        ## slhafile = installDirectory()+'oldFiles/andrePT4.slha'
        slhafile = '../inputFiles/slha/lightSquarks.slha'
        self.configureLogger()
        smstoplist = slhaDecomposer.decompose(slhafile, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        listofanalyses = smsAnalysisFactory.load()
        for analysis in listofanalyses:
            self.checkAnalysis(analysis,smstoplist)

if __name__ == "__main__":
    unittest.main()
