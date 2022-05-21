#!/usr/bin/env python3

"""
.. module:: testRuntime
   :synopsis: Tests the runtime methods
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""
import sys
sys.path.insert(0,"../")
from smodels.tools import runtime
import unittest

class TxTest(unittest.TestCase):
    from smodels.tools.smodelsLogging import logger

    def testNCPUs(self):
        ncpus = runtime.nCPUs()
        self.assertTrue ( ncpus >= 1 )

    def testDetermineNCPUs ( self ):
        from smodels.tools.modelTester import _determineNCPus
        ncpus = runtime.nCPUs()
        n = _determineNCPus ( 0, 1e6 )
        self.assertTrue ( n == ncpus )
        n = _determineNCPus ( -1e16, 1e6 )
        self.assertTrue ( n == 1 )

if __name__ == "__main__":
    unittest.main()
    # a=TxTest("testNCPUs") 
    # a.testNCPUs()
    # a.debug()
    # a.run()
