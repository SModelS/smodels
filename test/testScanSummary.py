#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import compareSummaries, runMain

class ScanSummaryTest(unittest.TestCase):

    def testPythonSummary( self ):
        out = "./unitTestOutput"
        dirname = "./testFiles/slha/"
        runMain(dirname,inifile = "testParameters.ini")
        outSummary = os.path.join(out,'summary.txt')
        outDefault = 'summary_scan_default.txt'

        self.assertTrue(compareSummaries(outSummary,outDefault,allowedDiff=0.05))


    def testSLHASummary( self ):
        out = "./unitTestOutput"
        dirname = "./testFiles/slha/"
        runMain(dirname,inifile = "slhaOnly.ini",suppressStdout=False,
                             development=True)
        outSummary = os.path.join(out,'summary.txt')
        outDefault = 'summary_scan_default.txt'
        self.assertTrue(compareSummaries(outSummary,outDefault,allowedDiff=0.05))

    def testSummarySummary( self ):
        out = "./unitTestOutput"
        for i in os.listdir( out ):
            if i[-8:]==".smodels":
                os.unlink(os.path.join(out, i))
        dirname = "./testFiles/slha/"
        runMain(dirname,inifile = "summaryOnly.ini")
        outSummary = os.path.join(out,'summary.txt')
        outDefault = 'summary_scan_default.txt'
        ret = compareSummaries(outSummary,outDefault,allowedDiff=0.05)
        if not ret:
            print ( f"ERROR: files {outSummary} and {outDefault} differ!" )
        self.assertTrue(ret )

if __name__ == "__main__":
    unittest.main()
