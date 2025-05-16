#!/usr/bin/env python3

"""
.. module:: testPrinter
   :synopsis: Tests printer facilities with runSModelS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import os
import subprocess
sys.path.insert(0, "../")
import unittest
from xml.etree import ElementTree
from unitTestHelpers import equalObjs, Summary, runMain, importModule, sortXML,compareXML,compareSLHA
from smodels.base.smodelsLogging import logger
inf = float("inf")



class RunPrinterTest(unittest.TestCase):
    definingRun = False  # meant only to adapt to changes in output format
    # use with super great care!!

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [f, f.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)

    def testPrintersV2(self):

        slhafile = "./testFiles/slha/gluino_squarks.slha"
        out = runMain(slhafile)
        outputfile = out.replace('.py', '.smodels')

        defaultfile = "gluino_squarks_default.smodels"
        # Test summary output
        output = Summary(outputfile, allowedRelDiff=0.05)
        default = Summary(defaultfile, allowedRelDiff=0.05)
        if default != output:
            logger.error ( f"{outputfile} and {defaultfile} differ!" )

        self.assertEqual(default, output)
        self.removeOutputs(outputfile)

        smodelsOutput = importModule(out)
        # Test python output
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields, where="top",
                           fname="./unitTestOutput/printer_output.py")
        if default != output:
            logger.error ( f"{outputfile} and gluino_squarks_default.py differ!" )

        self.assertTrue(equals)
        self.removeOutputs(out)
        self.removeOutputs('./debug.log')  

        outputfile = out.replace('.py', '.xml')
        defFile = "default_output.xml"
        # Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        comp = compareXML(xmlDefault, xmlNew,
                      allowedRelDiff=0.05,
                      ignore=['input_file', 'smodels_version', 'database_version', 'ncpus'])
        if not comp:
            logger.error ( f"{outputfile} and {defFile} differ!" )
        self.assertTrue(comp)
        self.removeOutputs(outputfile)
        self.removeOutputs('./debug.log')

        outputfile = out.replace('.py', '.smodelsslha')

        slhaDefaultFile = "./gluino_squarks_default.slha.smodelsslha"
        comp = compareSLHA(slhaDefaultFile, outputfile,
                           allowedRelDiff=0.05)
        if not comp:
            logger.error ( f"{outputfile} and {slhaDefaultFile} differ!" )
        self.assertTrue(comp )
        if comp:
            self.removeOutputs(outputfile)

    def testPythonPrinterSimpleV2(self):

        slhafile = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(slhafile, inifile='testParameters_exp.ini')
        smodelsOutput = importModule(outputfile)

        if self.definingRun:
            from smodels.base.smodelsLogging import logger
            logger.error("This is a definition run! Know what youre doing!")
            default = "simplyGluino_default.py"
            outputfile = './unitTestOutput/printer_output_simple.py'
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % (outputfile, default)
            subprocess.getoutput(cmd)
        from simplyGluino_default_extended import smodelsOutputDefault

        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version', 'model', 'checkinput',
                        'doinvisible', 'docompress', 'computestatistics', 'testcoverage']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields)
        from smodels.base.smodelsLogging import logger

        if equals == False:
            logger.error ( f"{outputfile} and simplyGluino_default_extended.py differ!" )
        self.assertTrue(equals)
        self.removeOutputs(outputfile)

        outputfile = outputfile.replace('.py', '.xml')

        defFile = "default_outputSimplyGluino.xml"

        # Test xml output
        xmlDefault = ElementTree.parse(defFile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        ret = compareXML(xmlDefault, xmlNew, allowedRelDiff=0.05, 
               ignore=['input_file', 'smodels_version', 'database_version', 'ncpus' ] )
        if not ret:
            logger.error ( f"difference between {defFile} and {outputfile}" )

        self.assertTrue( ret)

        if ret:
            self.removeOutputs(outputfile)
            self.removeOutputs('./debug.log')


    def testPrinters(self):

        slhafile = "./testFiles/slha/lightEWinos.slha"
        out = runMain(slhafile,inifile="testPrinters_parameters.ini",
                suppressStdout = False )

        # Check Summary output
        outputfile = out.replace('.py', '.smodels')
        defaultfile = "lightEWinos_default.smodels"
        output = Summary(outputfile, allowedRelDiff=0.05)
        default = Summary(defaultfile, allowedRelDiff=0.05)
        if output != default:
            logger.error ( f"{out} differs from {defaultfile}" )
        else:
            self.assertEqual(default, output)
            self.removeOutputs(outputfile)

        # Check Python output
        smodelsOutput = importModule(out)
        from lightEWinos_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields, where="top")
        if not equals:
            logger.error ( f"{out} differs from lightEWinos_default.py" )
            
        self.assertTrue(equals)
        if equals:
            self.removeOutputs(out)
            self.removeOutputs('./debug.log')       


        # Check XML output:
        outputfile = out.replace('.py', '.xml')
        defaultfile = "lightEWinos_default.xml"
        xmlDefault = ElementTree.parse(defaultfile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        equals = compareXML(xmlDefault, xmlNew, allowedRelDiff=0.05,
                            ignore=['input_file', 'smodels_version', 'database_version', 'ncpus'])
        if not equals:
            logger.error ( f"{outputfile}!={defaultfile}" )

        self.assertTrue( equals )
        if equals:
            self.removeOutputs(outputfile)
            self.removeOutputs('./debug.log')
        

        # Check SLHA output:
        outputfile = out.replace('.py', '.smodelsslha')
        slhaDefaultFile = "./lightEWinos_default.smodelsslha"
        comp = compareSLHA(slhaDefaultFile, outputfile,
                           allowedRelDiff=0.05)
        if not comp:
            logger.error ( f"{outputfile} and {slhaDefaultFile} differ!" )
        self.assertTrue(comp)
        if comp:
            self.removeOutputs(outputfile)

    
    def testPythonPrinterNodesMap(self):

        slhafile = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(slhafile, inifile='testPrinters_parameters_nodeMap.ini')
        smodelsOutput = importModule(outputfile)
        from simplyGluino_default_nodesMap import smodelsOutputDefault

        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version', 'model', 'checkinput',
                        'doinvisible', 'docompress', 'computestatistics', 'testcoverage']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        smodelsOutput['ExptRes'] = sorted(smodelsOutput['ExptRes'],
                                          key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields)

        self.assertTrue(equals)
        self.removeOutputs(outputfile)


if __name__ == "__main__":
    unittest.main()
