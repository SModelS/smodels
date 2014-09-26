#!/usr/bin/env python

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

from smodels.tools import toolBox
import unittest
import logging

class TestToolBox(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def testToolBox(self):
        self.logger.info ( "ToolBox" )
        box = toolBox.ToolBox()
        ok=box.checkInstallation ( colors=False, printit=False )
        self.assertEqual ( ok, True )

if __name__ == "__main__":
    unittest.main()
