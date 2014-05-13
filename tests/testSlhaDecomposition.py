#!/usr/bin/env python

"""
.. module:: testSlhaDecomposition
   :synopsis: Checks slha decomposition, alongside with compression
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import setPath
from smodels.theory import slhaDecomposer
import unittest
import logging
#import logging.config

class SlhaDecompositionTest(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def test(self):
        self.logger.info ( "test decomposition, no compression" )
        """ test the decomposition with no compression """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, 0.1, False, False )
        for topo in topos:
            print topo
            for element in topo.elementList:
                print element,element.mother

if __name__ == "__main__":
    unittest.main()
