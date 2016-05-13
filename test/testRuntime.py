#!/usr/bin/env python

"""
.. module:: testRuntime
   :synopsis: Tests the runtime methods
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""
import sys
sys.path.insert(0,"../")
from smodels.tools import runtime
import unittest
import logging

class TxTest(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def testNCPUs(self):
        ncpus = runtime.nCPUs()
        self.assertTrue ( ncpus >= 1 )

if __name__ == "__main__":
    unittest.main()
