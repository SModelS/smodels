#!/usr/bin/env python

"""
.. module:: testSlhaDecomposition
   :synopsis: Checks slha decomposition, alongside with compression
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import setPath
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV, fb
import unittest
import logging
#import logging.config

class SlhaDecompositionTest(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def test(self):
        self.logger.info ( "test decomposition, no compression" )
        """ test the decomposition with no compression """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, False, 5.*GeV )
        self.assertEqual ( len(topos), 2 )
        e1,e2=len(topos[0].elementList),len(topos[1].elementList)
        if e1>e2:
            e1,e2=e2,e1
        self.assertEqual ( e1, 12 )
        self.assertEqual ( e2, 120 )

if __name__ == "__main__":
    unittest.main()
