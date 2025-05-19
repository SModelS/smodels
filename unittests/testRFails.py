#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import warnings
import sys
import os
sys.path.insert(0, "../")
import unittest
from xml.etree import ElementTree
from unitTestHelpers import equalObjs, runMain, importModule, Summary, sortXML,compareXML, compareSLHA, compareScanSummary
from smodels.base.smodelsLogging import logger
from smodels.experiment.databaseObj import Database

class RFailsTest(unittest.TestCase):
    definingRun = False  # meant only to adapt to changes in output format
    #  use with super great care!!


    def removeOutputs(self, f):
        """ remove cruft outputfiles """

        f = os.path.splitext(f)[0]
        extList = ['py','pyc','smodels','smodelsslha','xml']
        for ext in extList:
            fname= f+'.'+ext
            if os.path.exists(fname):
                os.remove(fname)

    def testRFailPrinters(self):
        warnings.filterwarnings( action='ignore', category=DeprecationWarning )
        filename = "./testFiles/rfails/mdm_r_fails.slha"
        from smodels.statistics import pyhfInterface

        # Reduce the number of attempts 
        pyhfInterface.nattempts_max = 5
        database = Database('official')
        out = runMain(filename,inifile='testParameters_rfails.ini',
                             overridedatabase = database)

        # Check Summary output
        outputfile = out.replace('.py', '.smodels')
        defaultfile = "mdm_r_fails_default.smodels"
        output = Summary(outputfile, allowedRelDiff=0.05)
        default = Summary(defaultfile, allowedRelDiff=0.05)
        if output != default:
            logger.error ( f"{out} differs from {defaultfile}" )
        else:
            self.assertEqual(default, output)
            
        # Check Python output
        smodelsOutput = importModule(out)
        from mdm_r_fails_default import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version']
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.05,
                           ignore=ignoreFields, where="top")
        if not equals:
            logger.error ( f"{out} differs from mdm_r_fails_default.py" )
            
        self.assertTrue(equals)


        # Check XML output:
        outputfile = out.replace('.py', '.xml')
        defaultfile = "mdm_r_fails_default.xml"
        xmlDefault = ElementTree.parse(defaultfile).getroot()
        xmlNew = ElementTree.parse(outputfile).getroot()
        sortXML(xmlDefault)
        sortXML(xmlNew)
        equals = compareXML(xmlDefault, xmlNew, allowedRelDiff=0.05,
                            ignore=['input_file', 'smodels_version', 'database_version', 'ncpus'])
        if not equals:
            logger.error ( f"{outputfile}!={defaultfile}" )

        self.assertTrue( equals )
        

        # Check SLHA output:
        outputfile = out.replace('.py', '.smodelsslha')
        slhaDefaultFile = "./mdm_r_fails_default.smodelsslha"
        equals = compareSLHA(slhaDefaultFile, outputfile,
                           allowedRelDiff=0.05)
        if not equals:
            logger.error ( f"{outputfile} and {slhaDefaultFile} differ!" )
        self.assertTrue(equals)
        if equals:
            self.removeOutputs(outputfile)

    def testRFailSummary(self):
        warnings.filterwarnings( action='ignore', category=DeprecationWarning )
        warnings.filterwarnings( action='ignore', category=RuntimeWarning )
        out = "./unitTestOutput"
        dirname = "./testFiles/rfails/"
        from smodels.statistics import pyhfInterface

        # Reduce the number of attempts 
        pyhfInterface.nattempts_max = 10
        database = Database('official')
        runMain(dirname,inifile='testParameters_rfails.ini',
                             overridedatabase = database,suppressStdout=True )
        outSummary = os.path.join(out, 'summary.txt')
        outDefault = 'r_fails_summary_default.txt'
        comp = compareScanSummary(outSummary, outDefault, allowedRelDiff=0.05)
        if not comp:
            print ( f"ERROR: {outSummary}!={outDefault}" )

        self.assertTrue(comp)

        if comp:
            for f in os.listdir(dirname):
                self.removeOutputs(os.path.join(out,os.path.basename(f)))

            if os.path.isfile(outSummary):
                os.remove(outSummary)

if __name__ == "__main__":
    unittest.main()
