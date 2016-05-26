#!/usr/bin/env python

"""
.. module:: testPrinter
   :synopsis: Tests printer facilities with runSModelS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory as idir
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import fb, GeV, TeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools import printer, ioObjects
from smodels.tools import coverage
from smodels.tools import summaryReader
from xml.etree import ElementTree
import logging
from databaseLoader import database
from os.path import join

class RunPrinterTest(unittest.TestCase):

    def __init__ ( self, *args, **kwargs):
        super(RunPrinterTest, self).__init__(*args, **kwargs)
        self.logger = logging.getLogger(__name__)
        self.stdoutPrinter = printer.TxTPrinter(output = 'file', 
                filename = './unitTestOutput/sms_output.txt')
        self.summaryPrinter = printer.SummaryPrinter(output = 'file', 
                filename = './unitTestOutput/summary_print.txt')
        self.pythonPrinter = printer.PyPrinter(output = 'file', 
                filename = './unitTestOutput/sms_output.py')
        self.xmlPrinter = printer.XmlPrinter(output = 'file', 
                filename = './unitTestOutput/sms_output.xml')
        self.printerList = printer.MPrinter( self.stdoutPrinter,self.summaryPrinter,
                                             self.pythonPrinter,self.xmlPrinter)
        #Set the address of the database folder
        self.slhafile = join ( idir(), "inputFiles/slha/gluino_squarks.slha" )
        self.runMain()


    def runMain(self):
        """
        Main program. Displays basic use case.
    
        """

        #Set main options for decomposition:
        sigmacut = 0.03 * fb
        mingap = 5. * GeV
    
        """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer
            for LHE input) """
        smstoplist = slhaDecomposer.decompose(self.slhafile, sigmacut, 
                         doCompress=True, doInvisible=True, minmassgap=mingap )
    
        #Add the decomposition result to the printers
        self.printerList.addObj(smstoplist)
    
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
        self.printerList.addObj(results,objOutputLevel=2)
        
        #Add coverage information:
        sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
        coverageInfo = coverage.Uncovered(smstoplist, False, False, sqrts)
        self.printerList.addObj(coverageInfo.missingTopos)
        
        #Add additional information:
        databaseVersion = database.databaseVersion
        outputStatus = ioObjects.OutputStatus([1,None], self.slhafile,
                                               {'sigmacut' : sigmacut.asNumber(fb), 
                                                'minmassgap' : mingap.asNumber(GeV),
                                                'maxcond': maxcond },
                                              databaseVersion)
        outputStatus.status = 1
        self.printerList.addObj(outputStatus)
        self.printerList.close()

    def describeDifferences ( self, A, B, defFile, testFile ):
        msg = "Dictionaries in %s and %s are different!" % ( defFile, testFile )
        Akeys=A.keys()
        Bkeys=B.keys()
        Akeys.sort()
        Bkeys.sort()
        if Akeys != Bkeys:
            self.logger.error ( "difference in keys:" )
            self.logger.error ( Akeys )
            self.logger.error ( Bkeys )
            raise AssertionError ( msg )
        for k in Akeys:
            Alines = A[k]
            Blines = B[k]
            if type ( Alines ) != type([]):
                if Alines != Blines:
                    self.logger.error ( "k=%s:default=%s, unittest=%s" % ( k, Alines, Blines ) )
                    raise AssertionError ( msg )
            else:
                if len(Alines)!=len(Blines):
                    self.logger.error ( "k=%s: # lines: %d != %d" % ( k, len(Alines), len(Blines) )  )
                    raise AssertionError ( msg )
                else:
                    raiseError=False
                    for Aline,Bline in zip(Alines, Blines):
                        if type(Aline)==type({}):
                            for kk in Aline.keys():
                                Akv = Aline[kk]
                                Bkv = Bline[kk]
                                if type ( Akv ) == float:
                                    if abs ( Akv - Bkv ) > 10**-7:
                                        self.logger.error ( "k=%s: kk=%s: %.13f != %.13f" % ( k, kk, Akv, Bkv ) )
                                        raiseError=True
                                else:
                                    if Akv != Bkv:
                                        self.logger.error ( "k=%s: kk=%s: %s != %s" % ( k, kk, Akv, Bkv ) )
                                        raiseError=True
                            continue
                        if Aline!=Bline:
                            self.logger.error ( Aline )
                            self.logger.error ( Bline )
                            raiseError=True
                    if raiseError:
                        raise AssertionError ( msg )

    def testTextPrinter(self):
        outputfile = join ( idir(), "test/unitTestOutput/summary_print.txt" )
        samplefile = join ( idir(), "test/summary_default.txt" )
        #Test summary output
        output = summaryReader.Summary( outputfile )
        sample = summaryReader.Summary( samplefile )
        try:
            self.assertEquals(sample, output)
        except AssertionError,e:
            msg = "%s != %s" % ( sample, output ) 
            raise AssertionError ( msg )


    def testPythonPrinter(self):
        #Test python output
        from default_output import smodelsOutput as defaultOut
        from unitTestOutput.sms_output import smodelsOutput
        inputFile = smodelsOutput['input file']
        inputFile = inputFile[inputFile.rfind('/'):]
        smodelsOutput['input file'] = inputFile
        defaultFile = defaultOut['input file']
        defaultFile = defaultFile[defaultFile.rfind('/'):]
        defaultOut['input file'] = defaultFile
        self.describeDifferences ( smodelsOutput, defaultOut, "./default_output.py", "./unitTestOutput/sms_output.py" )

    def testXmlPrinter(self):
        defFile = join ( idir(), "test/default_output.xml" )
        outFile = join ( idir(), "test/unitTestOutput/sms_output.xml" )
        
        #Test xml output
        xmlDefault = ElementTree.parse( defFile ).getroot()
        xmlNew = ElementTree.parse( outFile ).getroot()
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
        try:
            compareXML(xmlDefault,xmlNew)
        except AssertionError,e:
            msg = "%s != %s" % ( defFile, outFile )
            raise AssertionError ( msg )


if __name__ == "__main__":
    unittest.main()
