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
from smodels.base import runtime

class SpeyTest(unittest.TestCase):
    def testSwitch(self):
        """ see that we can turn on spey mode """
        computer = getStatsComputerModule()
        self.assertTrue( type(computer) == type(StatsComputer) )
        runtime._experimental["spey"]=True
        computer = getStatsComputerModule()
        self.assertTrue( type(computer) == type(SpeyComputer) )
        ## important! need to set back
        runtime._experimental["spey"]=False

if __name__ == "__main__":
    unittest.main()
