#!/usr/bin/env python3

"""
.. module:: testSpey
   :synopsis: Test the spey hidden feature

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys

sys.path.insert(0, "../")
import unittest

# from smodels.tools import statistics
from smodels.statistics.statsTools import getStatsComputerModule, StatsComputer
from unitTestHelpers import equalObjs, runMain, importModule, removeCruftOutputs
from smodels.base import runtime

class SpeyTest(unittest.TestCase):
    def testSwitch(self):
        from smodels.statistics.speyTools import SpeyComputer
        """ see that we can turn on spey mode """
        computer = getStatsComputerModule()
        self.assertTrue( type(computer) == type(StatsComputer) )
        runtime._experimental["spey"]=True
        computer = getStatsComputerModule()
        self.assertTrue( type(computer) == type(SpeyComputer) )
        ## important! need to set back
        runtime._experimental["spey"]=False

    def testIniFile(self):
        """ see that we can turn on spey mode """
        filename = "./testFiles/slha/gluino_squarks.slha"
        inifile = "testParameters_spey.ini"
        outputfile = runMain(filename, inifile = inifile, suppressStdout=False )
        smodelsOutput = importModule(outputfile)
        from default_with_spey import smodelsOutputDefault
        runtime._experimental["spey"]=False
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'Element',
                        'database version', 'model']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                           key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, 
                           allowedRelDiff=0.02,
                           ignore=ignoreFields, fname=outputfile)
        if not equals:
            p = outputfile.find("unitTestOutput")
            fname = outputfile
            if p > 0:
                fname = fname[p:]
            print(f"[testRunSModelS] {fname} != gluino_squarks_default.py")
        self.assertTrue(equals)
        removeCruftOutputs(outputfile)

if __name__ == "__main__":
    unittest.main()
