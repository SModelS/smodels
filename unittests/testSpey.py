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
from smodels.statistics.statsTools import getCompRetrieverModule, \
         StatsComputer, CompRetriever
from unitTestHelpers import equalObjs, runMain, importModule, removeCruftOutputs
from smodels.base import runtime
import warnings

class SpeyTest(unittest.TestCase):
    def testSwitch(self):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message=r".*jsonschema\.RefResolver is deprecated.*"
            )
            from smodels.statistics.speyTools import SpeyRetriever
            """ see that we can turn on spey mode """
            computer = getCompRetrieverModule()
            self.assertTrue( type(computer) == type(CompRetriever) )
            runtime._experimental["spey"]=True
            computer = getCompRetrieverModule()
            self.assertTrue( type(computer) == type(SpeyRetriever) )
            ## important! need to set back
            runtime._experimental["spey"]=False

    def testIniFile(self):
        """ see that we can turn on spey mode """
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message=r".*jsonschema\.RefResolver is deprecated.*"
            )
            filename = "./testFiles/slha/gluino_squarks.slha"
            inifile = "testParameters_spey.ini"
            from databaseLoader import database
            outputfile = runMain( filename, inifile = inifile,
                                  overridedatabase = database,
                                  suppressStdout=True )
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
                print(f"[testRunSModelS] {fname} != default_with_spey.py")
            self.assertTrue(equals)
            removeCruftOutputs(outputfile)

if __name__ == "__main__":
    unittest.main()
