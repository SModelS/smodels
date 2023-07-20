#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0, "../")
import os
import unittest
from unitTestHelpers import compareScanSummary, runMain
from smodels.base.smodelsLogging import logger, setLogLevel, getLogLevel
setLogLevel("error")

class ScanSummaryTest(unittest.TestCase):

    # Note that in v3 in get results for the idm_example.slha file, because
    # the heavy higgses are part of the MSSM model file and we no longer
    # have the Z2odd requirement.

    def removeOutputs(self,f):
        extList = ['py','pyc','smodels','smodelsslha','xml']
        for ext in extList:
            fname= f+'.'+ext
            if os.path.exists(fname):
                os.remove(fname)

    def testPythonSummary(self):
        out = "./unitTestOutput"
        dirname = "./testFiles/slha/"
        runMain(dirname, inifile="testParameters_v2.ini")
        outSummary = os.path.join(out, 'summary.txt')
        outDefault = 'summary_scan_default.txt'
        comp = compareScanSummary(outSummary, outDefault, allowedRelDiff=0.05)
        if not comp:
            print ( f"ERROR: {outSummary}!={outDefault}" )

        self.assertTrue(comp)

        for f in os.listdir(dirname):
            self.removeOutputs(os.path.join(out,os.path.basename(f)))

        if os.path.isfile(outSummary):
            os.remove(outSummary)

    def testSLHASummary(self):
        out = "./unitTestOutput"
        dirname = "./testFiles/slha/"
        runMain(dirname, inifile="slhaOnly.ini", suppressStdout=False,
                development=True)
        outSummary = os.path.join(out, 'summary.txt')
        outDefault = 'summary_scan_default.txt'
        comp = compareScanSummary(outSummary, outDefault, allowedRelDiff=0.05)
        if not comp:
            print ( f"ERROR: {outSummary}!={outDefault}" )
        self.assertTrue(comp)

        for f in os.listdir(dirname):
            self.removeOutputs(os.path.join(out,os.path.basename(f)))

        if os.path.isfile(outSummary):
            os.remove(outSummary)

    def testSummarySummary(self):
        out = "./unitTestOutput"
        for i in os.listdir(out):
            if i[-8:] == ".smodels":
                os.unlink(os.path.join(out, i))
        dirname = "./testFiles/slha/"
        runMain(dirname, inifile="summaryOnly.ini")
        outSummary = os.path.join(out, 'summary.txt')
        outDefault = 'summary_scan_default.txt'
        ret = compareScanSummary(outSummary, outDefault, allowedRelDiff=0.05)
        if not ret:
            print(f"ERROR: files {outSummary} and {outDefault} differ!")
        self.assertTrue(ret)


        for f in os.listdir(dirname):
            self.removeOutputs(os.path.join(out,os.path.basename(f)))

        if os.path.isfile(outSummary):
            os.remove(outSummary)



if __name__ == "__main__":
    unittest.main()
