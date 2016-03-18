#!/usr/bin/env python

"""
.. module:: testIntegration
   :synopsis: Integration test, tests a simple but complete use case.
              Uses the database in smodels/test/database.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.physicsUnits import fb, GeV, pb
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
        # print "checking expresult",expresult
        from smodels.theory.theoryPrediction import theoryPredictionsFor
        theorypredictions = theoryPredictionsFor(expresult, smstoplist)
        defpreds=self.predictions()
        #print "ana",expresult,theorypredictions
        if not theorypredictions:
            print "no theory predictions for",expresult,"??"
            import sys
            sys.exit(-1)
        #print(">>>> Ana %s" % expresult)
        for pred in theorypredictions:
            #print ( "Pred ana",str(expresult.info.getInfo('id')))
            #print ( "Pred txname",pred.txname.getInfo("txname"))
            m0=str ( int ( pred.mass[0][0]/GeV )  )
            #print ( "Pred mass0",m0 )
            #print ( "Pred mass0",int(pred.mass[0]/GeV))
            #print ( "Pred value",pred.value )
            #print ( "Pred value type",type(pred.value) )
            #print ( "Pred value dict",pred.value.__dict__ )
            w=pred.value.getDictionary()[(None,None)]
            predval=w['8 TeV (NLL)']
            #print w,"%10f" % (predval/fb)
            #print "predval=",predval.asNumber(fb)
            defpredval=defpreds[expresult.getValuesFor('id')[0]+":"+pred.txnames[0].getInfo("txName") ]
            #print ( "Pred default value",defpredval )
            self.assertAlmostEqual ( predval / fb, defpredval / fb )

            # print ( "Pred condition",pred.conditions)
            ## print ( "Pred mass",pred.mass)

    def testIntegration(self):
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import fb, GeV
        from smodels.theory import slhaDecomposer
        from smodels.experiment.databaseObj import Database
        #from smodels.experiment import smsAnalysisFactory, smsHelpers
        #smsHelpers.base = installDirectory() + 'test/database/'
        # smsHelpers.runs = [ "2012" ]
        ## slhafile = installDirectory()+'oldFiles/andrePT4.slha'
        #slhafile = '../inputFiles/slha/lightSquarks.slha'
        slhafile = '../inputFiles/slha/simplyGluino.slha'
        ## slhafile = '../inputFiles/slha/compression.slha'
        self.configureLogger()
        smstoplist = slhaDecomposer.decompose(slhafile, .1*fb, doCompress=True,
                doInvisible=True, minmassgap=5.*GeV)
        database = Database ( "./database/" )
        listofanalyses = database.getExpResults( analysisIDs= [ "ATLAS-SUSY-2013-02" ], txnames = [ "T1" ] )
        ## print "analyses=",listofanalyses
        if type(listofanalyses) != list:
            listofanalyses= [ listofanalyses] 
        for analysis in listofanalyses:
            self.checkAnalysis(analysis,smstoplist)

if __name__ == "__main__":
    unittest.main()
