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
from smodels.tools.physicsUnits import fb, GeV, TeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.tools import printer, ioObjects
from smodels.tools import coverage
from smodels.tools import summaryReader
from xml.etree import ElementTree

stdoutPrinter = printer.TxTPrinter(output = 'file', filename = './unitTestOutput/sms_output.txt')
summaryPrinter = printer.SummaryPrinter(output = 'file', filename = './unitTestOutput/summary_print.txt')
pythonPrinter = printer.PyPrinter(output = 'file', filename = './unitTestOutput/sms_output.py')
xmlPrinter = printer.XmlPrinter(output = 'file', filename = './unitTestOutput/sms_output.xml')
printerList = printer.MPrinter(stdoutPrinter,summaryPrinter,pythonPrinter,xmlPrinter)
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
        
        #Add coverage information:
        sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
        coverageInfo = coverage.Uncovered(smstoplist, False, False, sqrts)
        printerList.addObj(coverageInfo.missingTopos)
        
        #Add additional information:
        databaseVersion = database.databaseVersion
        outputStatus = ioObjects.OutputStatus([1,None], slhafile,
                                               {'sigmacut' : sigmacut.asNumber(fb), 
                                                'minmassgap' : mingap.asNumber(GeV),
                                                'maxcond': maxcond },
                                              databaseVersion)
        outputStatus.status = 1
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
        try:
            self.assertEquals(sample, output)
        except AssertionError,e:
            raise AssertionError ( "%s != %s" % ( sample, output ) )
        
        #Test python output
        from default_output import smodelsOutput as defaultOut
        from unitTestOutput.sms_output import smodelsOutput
        inputFile = smodelsOutput['input file']
        inputFile = inputFile[inputFile.rfind('/'):]
        smodelsOutput['input file'] = inputFile
        defaultFile = defaultOut['input file']
        defaultFile = defaultFile[defaultFile.rfind('/'):]
        defaultOut['input file'] = defaultFile
        try:
            self.assertEquals(smodelsOutput,defaultOut)
        except AssertionError,e:
            msg = "%s != %s" % ( smodelsOutput, defaultOut )
            raise AssertionError ( msg )
        
        #Test xml output
        xmlDefault = ElementTree.parse("%s/test/default_output.xml" %installDirectory()).getroot()
        xmlNew = ElementTree.parse("%s/test/unitTestOutput/sms_output.xml" %installDirectory()).getroot()
        def compareXML(xmldefault,xmlnew):
            self.assertEqual(len(xmldefault),len(xmlnew))
            for i,el in enumerate(xmldefault):
                newel = xmlnew[i]
                self.assertEqual(len(el),len(newel))                
                if len(el) == 0:
                    if el.tag == 'input_file': continue
                    self.assertEqual(el.text,newel.text)
                    self.assertEqual(el.tag,newel.tag)
                else:                    
                    compareXML(el,newel)
        compareXML(xmlDefault,xmlNew)
        

if __name__ == "__main__":
    unittest.main()
