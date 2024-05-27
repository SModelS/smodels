#!/usr/bin/env python3

"""
.. module:: testReportAll
   :synopsis: Tests the reporting of all datasets

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.append('../')
import os
import unittest
from unitTestHelpers import equalObjs, runMain, importModule, sortSModelSOutput


class TestReportAll(unittest.TestCase):
    definingRun = False  # meant only to adapt to changes in output format
    # use with super great care!!

    def removeOutputs(self, f):
        """ remove cruft outputfiles """

        f = os.path.splitext(f)[0]
        extList = ['py','pyc','smodels','smodelsslha','xml']
        for ext in extList:
            fname= f+'.'+ext
            if os.path.exists(fname):
                os.remove(fname)

    def testReportAll(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        inifile = "testReportAll.ini"
        outputfile = runMain(filename, inifile=inifile, suppressStdout=True)
        smodelsOutput = importModule(outputfile)
        from default_reportAll import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus',
                        'database version']
        smodelsOutput = sortSModelSOutput ( smodelsOutput )
        smodelsOutputDefault = sortSModelSOutput ( smodelsOutputDefault )
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.03,
                           ignore=ignoreFields, fname=outputfile,
                           fname2="default_reportAll.py")
        self.assertTrue(equals)
        if equals:
            self.removeOutputs(outputfile)


if __name__ == "__main__":
    unittest.main()
