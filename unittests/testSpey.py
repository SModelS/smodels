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
from smodels.statistics.speyTools import SpeyComputer
from unitTestHelpers import equalObjs, runMain
from smodels.base import runtime

class SpeyTest(unittest.TestCase):
    def restSwitch(self):
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
        print ( f"output {outputfile}" )
        cmd = f"cp {outputfile} bla.py"
        import subprocess
        subprocess.getoutput ( cmd )

if __name__ == "__main__":
    unittest.main()
