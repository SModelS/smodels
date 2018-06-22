#!/usr/bin/env python3

"""
.. module:: testPrinter
   :synopsis: Tests printer facilities with runSModelS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory as idir
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.smodelsLogging import logger
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools import printer, ioObjects
from smodels.tools import coverage
import summaryReader
from xml.etree import ElementTree
from databaseLoader import database
import unum

def equalObjs(obj1,obj2,allowedDiff,ignore=[], where=None ):
    """
    Compare two objects.
    The numerical values are compared up to the precision defined by allowedDiff.
    
    :param obj1: First python object to be compared 
    :param obj2: Second python object to be compared
    :param allowedDiff: Allowed % difference between two numerical values 
    :param ignore: List of keys to be ignored
    :return: True/False    
    """      

    if type(obj1) in [ float, int ] and type ( obj2) in [ float, int ]:
        obj1,obj2=float(obj1),float(obj2)
    
    if type(obj1) != type(obj2):
        logger.info("Data types differ (%s,%s)" %(type(obj1),type(obj2)))
        return False
    
    if isinstance(obj1,unum.Unum):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff.asNumber() < allowedDiff
    elif isinstance(obj1,float):        
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff < allowedDiff
    elif isinstance(obj1,str):
        return obj1 == obj2
    elif isinstance(obj1,dict):    
        for key in obj1:
            if key in ignore: continue
            if not key in obj2:
                logger.warning("Key %s missing" %key)
                return False
            if not equalObjs(obj1[key],obj2[key],allowedDiff,ignore=ignore, where=key ):
                logger.warning('Objects differ in %s:\n   %s\n and\n   %s' %(where, str(obj1[key]),str(obj2[key])))
                return False
    elif isinstance(obj1,list):
        if len(obj1) != len(obj2):
            logger.warning("Objects differ in length:\n %s (len = %i) \n and\n %s (len = %i)" 
                            %(str(obj1),len(obj1),str(obj2),len(obj2)))
            return False
        for ival,val in enumerate(obj1):
            if not equalObjs(val,obj2[ival],allowedDiff, where = where ):
                logger.warning('Objects differ:\n   %s \n and\n   %s' %( str(val),str(obj2[ival])))
                return False
    else:
        return obj1 == obj2
            
    return True



def sortXML(xmltree):
    for el in xmltree:        
        sortXML(el)
    xmltree[:] = sorted(xmltree, key=lambda el: [el.tag,ElementTree.tostring(el)])

def compareXML(xmldefault,xmlnew,allowedDiff,ignore=[]):
    if len(xmldefault) != len(xmlnew):
        return False
    for i,el in enumerate(xmldefault):
        newel = xmlnew[i]
        if len(el) != len(newel):
            return False                
        if len(el) == 0:
            if el.tag in ignore: continue
            try:
                el.text = eval(el.text)
                newel.text = eval(newel.text)
            except (TypeError,NameError,SyntaxError):
                pass
            if isinstance(el.text,float) and isinstance(newel.text,float) \
                    and newel.text != el.text:
                diff = 2.*abs(el.text-newel.text)/abs(el.text+newel.text)
                if diff > allowedDiff:
                    return False
            else:
                if el.text != newel.text:
                    return False
            if el.tag != newel.tag:
                return False
        else:                    
            compareXML(el,newel,allowedDiff,ignore)

    return True



class RunPrinterTest(unittest.TestCase):

    def runPrinterMain(self, slhafile, mprinter,addTopList=False):
        """
        Main program. Displays basic use case.
    
        """

        #Set main options for decomposition:
        sigmacut = 0.03 * fb
        mingap = 5. * GeV
    
        """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer
            for LHE input) """
        smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, 
                         doCompress=True, doInvisible=True, minmassgap=mingap )
    
        #Add the decomposition result to the printers
        if addTopList:
            mprinter.addObj(smstoplist)
    
        listOfExpRes = database.getExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'])
        # Compute the theory predictions for each analysis
        allPredictions = []
        for expResult in listOfExpRes:
            predictions = theoryPredictionsFor(expResult, smstoplist)
            if not predictions:
                continue
            allPredictions += predictions._theoryPredictions
        
        maxcond = 0.2
        results = ioObjects.ResultList(allPredictions,maxcond)
        mprinter.addObj(results)
        
        #Add coverage information:
        coverageInfo = coverage.Uncovered(smstoplist)
        mprinter.addObj(coverageInfo)
        
        
        #Add additional information:
        databaseVersion = database.databaseVersion
        outputStatus = ioObjects.OutputStatus([1,'Input file ok'], slhafile,
                                               {'sigmacut' : sigmacut.asNumber(fb), 
                                                'minmassgap' : mingap.asNumber(GeV),
                                                'maxcond': maxcond },
                                              databaseVersion)
        outputStatus.status = 1
        mprinter.addObj(outputStatus)
        mprinter.flush()
        


    def testTextPrinter(self):
        outputfile = os.path.join( idir(), "test/unitTestOutput/printer_output.smodels")
        self.removeOutputs ( outputfile )
        mprinter = printer.MPrinter()
        mprinter.Printers['summary'] = printer.SummaryPrinter(output = 'file')        
        #Define options:
        mprinter.Printers['summary'].expandedSummary = True
         
        slhafile = os.path.join ( idir(), "inputFiles/slha/gluino_squarks.slha" )
        mprinter.setOutPutFiles('./unitTestOutput/printer_output',silent=True)
        self.runPrinterMain(slhafile,mprinter)
         
        samplefile = os.path.join( idir(), "test/gluino_squarks_default.txt")
        #Test summary output
        output = summaryReader.Summary(outputfile,allowedDiff=0.05)
        sample = summaryReader.Summary(samplefile,allowedDiff=0.05)
        try:
            self.assertEqual(sample, output)
        except AssertionError as e:
            msg = "%s != %s" %(sample, output) 
            raise AssertionError(msg)
        self.removeOutputs ( outputfile )
 
    def removeOutputs ( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists ( i ): os.remove ( i )
 
    def testPythonPrinter(self):
        self.removeOutputs ( './unitTestOutput/printer_output.py' )
        mprinter = printer.MPrinter()
        mprinter.Printers['python'] = printer.PyPrinter(output = 'file')
        #Define options:
        mprinter.Printers['python'].addElementList = False
               
        slhafile = os.path.join ( idir(), "inputFiles/slha/gluino_squarks.slha" )
        mprinter.setOutPutFiles('./unitTestOutput/printer_output',silent=True)
        self.runPrinterMain(slhafile,mprinter)
           
        from unitTestOutput.printer_output import smodelsOutput
        #Test python output
        from gluino_squarks_default import smodelsOutputDefault 
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version' ]
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'], 
                      key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                                       res['AnalysisID'],res['DataSetID']])
        equals = equalObjs( smodelsOutput,smodelsOutputDefault,allowedDiff=0.05,
                            ignore=ignoreFields, where = "top" )
        try:
            self.assertTrue(equals)
        except AssertionError as e:
            print ( "Error: %s, when comparing %s \nwith %s." % (e,"output.py","gluino_squarks_default.py" ) )
            raise AssertionError(e)
        self.removeOutputs ( './unitTestOutput/printer_output.py' )
 
    def testPythonPrinterSimple(self):
        self.removeOutputs ( './unitTestOutput/printer_output_simple.py' )
 
        mprinter = printer.MPrinter()
        mprinter.Printers['python'] = printer.PyPrinter(output = 'file')
        #Define options:
        mprinter.Printers['python'].addelementlist = True
         
        slhafile = os.path.join ( idir(), "inputFiles/slha/simplyGluino.slha" )
        mprinter.setOutPutFiles('./unitTestOutput/printer_output_simple',silent=True)
        self.runPrinterMain(slhafile,mprinter,addTopList=True)
        
        #Test python output
        from unitTestOutput.printer_output_simple import smodelsOutput
        from simplyGluino_default import smodelsOutputDefault    
         
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version' ]
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'], 
                       key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                                        res['AnalysisID'],res['DataSetID']])
        equals = equalObjs( smodelsOutput,smodelsOutputDefault,allowedDiff=0.05,
                            ignore=ignoreFields )
        try:
            self.assertTrue(equals)
        except AssertionError as e:
            print ( "Error: %s, when comparing %s \nwith %s." % (e,"outputSimple.py","simplyGluino_default.py" ) )
            raise AssertionError(e)
        self.removeOutputs ( './unitTestOutput/printer_output_simple.py' )
  
  
    def testXmlPrinter(self):
        self.removeOutputs ( './unitTestOutput/printer_output.xml' )  
  
        mprinter = printer.MPrinter()
        mprinter.Printers['xml'] = printer.XmlPrinter(output = 'file')
        #Define options:
        mprinter.Printers['xml'].addelementlist = False                    
  
        slhafile = os.path.join ( idir(), "inputFiles/slha/gluino_squarks.slha" )
        mprinter.setOutPutFiles('./unitTestOutput/printer_output',silent=True)
        self.runPrinterMain(slhafile,mprinter)                    
                      
        defFile = os.path.join ( idir(), "test/default_output.xml" )
        outFile = os.path.join ( idir(), "test/unitTestOutput/printer_output.xml" )
           
        #Test xml output
        xmlDefault = ElementTree.parse( defFile ).getroot()
        xmlNew = ElementTree.parse( outFile ).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        try:
            self.assertTrue(compareXML(xmlDefault,xmlNew,allowedDiff=0.05,ignore=['input_file','smodels_version', 'ncpus']))
        except AssertionError as e:
            msg = "%s != %s" %(defFile, outFile) + "\n" + str(e)            
            raise AssertionError(msg)
        self.removeOutputs ( './unitTestOutput/printer_output.xml' )  
 
    def testXmlPrinterSimple(self):
        self.removeOutputs ( './unitTestOutput/printer_output_simple.xml' )  
 
        mprinter = printer.MPrinter()
        mprinter.Printers['xml'] = printer.XmlPrinter(output = 'file')
        #Define options:
        mprinter.Printers['xml'].addelementlist = True   
 
        slhafile = os.path.join ( idir(), "inputFiles/slha/simplyGluino.slha" )
        mprinter.setOutPutFiles('./unitTestOutput/printer_output_simple',silent=True)
        self.runPrinterMain(slhafile,mprinter,addTopList=True)                    
                     
        defFile = os.path.join ( idir(), "test/default_outputSimplyGluino.xml" )
        outFile = os.path.join ( idir(), "test/unitTestOutput/printer_output_simple.xml" )
          
        #Test xml output
        xmlDefault = ElementTree.parse( defFile ).getroot()
        xmlNew = ElementTree.parse( outFile ).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        try:
            self.assertTrue(compareXML(xmlDefault,xmlNew,allowedDiff=0.05,ignore=['input_file','smodels_version', 'ncpus']))
        except AssertionError as e:
            msg = "%s != %s" %(defFile, outFile) + "\n" + str(e)            
            raise AssertionError(msg)
        self.removeOutputs ( './unitTestOutput/printer_output_simple.xml' )  


if __name__ == "__main__":
    unittest.main()
