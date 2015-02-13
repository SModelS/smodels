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
        self.assertEqual ( re, 18.0792727821 )

if __name__ == "__main__":
    unittest.main()
