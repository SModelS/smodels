#!/usr/bin/env python

"""
.. module:: testStatistics
   :synopsis: Tests the statistics module.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
from smodels.tools import statistics

class StatisticsTest(unittest.TestCase):
    def testCLInterval(self):
        re = statistics.computeCLInterval ( 100., 100., 1. )
        self.assertAlmostEqual ( re, 18.0792727821 )

    def testBayesianLimit(self):
        re = statistics.bayesianUpperLimit ( 100, 0., 100., 0. )
        self.assertAlmostEqual ( re, 21.4256127382 )

    def testUL(self):
        re = statistics.getUL ( 100, 100., 0. )
        print re

if __name__ == "__main__":
    unittest.main()
