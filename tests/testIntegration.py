#!/usr/bin/env python

"""
.. module:: testIntegration
   :synopsis: Integration test, tests a simple but complete use case.
              Uses the database in smodels/validation/database.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import setPath
from smodels.tools.physicsUnits import fb, GeV

class IntegrationTest(unittest.TestCase):
    def configureLogger(self):
        import logging.config
        logging.config.fileConfig( "/home/walten/git/smodels/tests/integration.conf",
               disable_existing_loggers=False )

    def predictions(self):
        return { 'SUS13006:TSlepSlep371': 2.42154060952*fb,
                 'SUS13006:TSlepSlep420': 0.246832973996*fb,
                 'SUS13006:TChiChipmSlepL371': 8.43114195825*fb,
                 'SUS12022:TChiChipmSlepL371': 8.43114195825*fb }

    def checkAnalysis(self,analysis,smstoplist):
        from smodels.theory.theoryPrediction import theoryPredictionFor
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        if not theorypredictions:
            return
        defpreds=self.predictions()
        ### print(">>>> Ana %s" % analysis.label)
        for pred in theorypredictions:
            # print ( "Pred ana",pred.analysis.label)
            m0=str ( int ( pred.mass[0][0]/GeV )  )
            # print ( "Pred mass0",m0 )
            # print ( "Pred mass0",int(pred.mass[0]/GeV))
            predval=pred.value.getDictionary()['8.0 [TeV]']
            defpredval=defpreds[pred.analysis.label+m0] 
            # print ( "Pred value",predval )
            # print ( "Pred default value",defpredval )
            self.assertAlmostEqual ( predval / fb, defpredval / fb )

            # print ( "Pred condition",pred.conditions)
            ## print ( "Pred mass",pred.mass)

    def testIntegration(self):
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import fb, GeV
        from smodels.theory import slhaDecomposer
        from smodels.experiment import smsAnalysisFactory, smsHelpers
        smsHelpers.base = installDirectory() + 'validation/database/'
        smsHelpers.runs = [ "2012" ]
        slhafile = installDirectory()+'inputFiles/slha/andrePT4.slha'
        self.configureLogger()
        smstoplist = slhaDecomposer.decompose(slhafile, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        listofanalyses = smsAnalysisFactory.load()
        for analysis in listofanalyses:
            self.checkAnalysis(analysis,smstoplist)

if __name__ == "__main__":
    unittest.main()
