#!/usr/bin/env python

"""
.. module:: testPrinter
   :synopsis: Tests printer facilities with runSModelS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.tools import printer, ioObjects
from smodels.tools import missingTopologies
from smodels.tools import summaryReader

stdoutPrinter = printer.TxTPrinter(output = 'file', filename = './unitTestOutput/sms_output.txt')
summaryPrinter = printer.SummaryPrinter(output = 'file', filename = './unitTestOutput/summary_print.txt')
pythonPrinter = printer.PyPrinter(output = 'file', filename = './unitTestOutput/sms_output.py')
printerList = printer.MPrinter(stdoutPrinter,summaryPrinter,pythonPrinter)
#Set the address of the database folder
database = Database("./database/")
slhafile = "%s/inputFiles/slha/gluino_squarks.slha" % \
                    (installDirectory() )



class RunPrinterTest(unittest.TestCase):
    def runMain(self):
        """
        Main program. Displays basic use case.
    
        """
    
        #Set main options for decomposition:
        sigmacut = 0.03 * fb
        mingap = 5. * GeV
    
        """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
        smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    
        #Add the decomposition result to the printers
        printerList.addObj(smstoplist)
    
        listOfExpRes = database.getExpResults()
        # Compute the theory predictions for each analysis
        allPredictions = []
        for expResult in listOfExpRes:
            predictions = theoryPredictionsFor(expResult, smstoplist)
            if not predictions:
                continue
            allPredictions += predictions._theoryPredictions
        
        maxcond = 0.2
        results = ioObjects.ResultList(allPredictions,maxcond)
        printerList.addObj(results,objOutputLevel=2)
        
        #Add missing topologies:
        sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
        missingtopos = missingTopologies.MissingTopoList(sqrts)
        missingtopos.findMissingTopos(smstoplist, listOfExpRes, mingap, True, doInvisible=True)        
        printerList.addObj(missingtopos)
        
        #Add additional information:
        databaseVersion = database.databaseVersion
        outputStatus = ioObjects.OutputStatus([1,None], slhafile,
                                               {'sigmacut (fb)' : sigmacut.asNumber(fb), 
                                                'mingap (GeV)' : mingap.asNumber(GeV)},
                                              databaseVersion)
        printerList.addObj(outputStatus)
        printerList.close()
        self.assertEqual(True, True)
        
    def testPrinters(self):
        self.runMain()
        
        #Test summary output
        output = summaryReader.Summary(
                "%s/test/unitTestOutput/summary_print.txt" %installDirectory())
        sample = summaryReader.Summary(
                "%s/test/summary_default.txt" %installDirectory())
        self.assertEquals(sample, output)
        
        #Test python output
        from default_output import smodelsOutput as defaultOut
        from unitTestOutput.sms_output import smodelsOutput
        inputFile = smodelsOutput['input file']
        inputFile = inputFile[inputFile.rfind('/'):]
        smodelsOutput['input file'] = inputFile
        defaultFile = defaultOut['input file']
        defaultFile = defaultFile[defaultFile.rfind('/'):]
        defaultOut['input file'] = defaultFile
        self.assertEquals(smodelsOutput,defaultOut)
        

if __name__ == "__main__":
    unittest.main()
