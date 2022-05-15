#!/usr/bin/env python3

"""
.. module:: testToolBox
   :synopsis: Tests the ToolBox thingie
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.tools import toolBox
import unittest

class TestToolBox(unittest.TestCase):
    from smodels.tools.smodelsLogging import logger

    def testToolBox(self):
        self.logger.info( "ToolBox" )
        box = toolBox.ToolBox()
        ok=box.checkInstallation( make=True, printit=False )
        self.assertTrue( ok )
 
if __name__ == "__main__":
    unittest.main()
