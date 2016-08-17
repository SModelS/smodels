#!/usr/bin/env python

"""
.. module:: testStatistics
   :synopsis: Tests the statistics module.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools import statistics
from smodels.tools.physicsUnits import fb
import scipy.stats
import math
import random

class StatisticsTest(unittest.TestCase):
    def testUpperLimit(self):
        re = statistics.upperLimit ( 100., 100., 0., 20./fb   )
        self.assertAlmostEqual ( re.asNumber ( fb ), 1.06, 1 )

if __name__ == "__main__":
    unittest.main()
