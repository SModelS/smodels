#!/usr/bin/env python

"""
.. module:: testToolBox
   :synopsis: Tests the ToolBox thingie
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.tools import toolBox
import unittest
import logging

class TestToolBox(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def testToolBox(self):
        self.logger.info ( "ToolBox" )
        box = toolBox.ToolBox()
        ok=box.checkInstallation ( colors=False, make=True, printit=False )
        self.assertTrue ( ok )

if __name__ == "__main__":
    unittest.main()
