#!/usr/bin/env python

"""
.. module:: testCompression
   :synopsis: Checks the compression algorithms
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import setPath
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV, fb
import unittest
import logging

class CompressionTest(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def testInvisiblePositive(self):
        """ test the invisible compression, a positive example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[0][1,1,0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[],[[nu],[ta+]]]":
                    continue
                # print element,"mother:",element.mother,element.compressionAlgorithm
                self.assertEqual ( str(element.mother),"[[[nu],[nu]],[[nu],[ta+]]]" )
                self.assertEqual ( str(element.compressionAlgorithm),"invisible" )

    def testInvisibleNegative(self):
        """ test the invisible compression, a negative example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[0][1,1,0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[],[[nu],[mu-]]]":
                    continue
                ## print element,"mother:",element.mother,element.compressionAlgorithm
                self.assertEqual ( str(element.mother),"None" )
                self.assertEqual ( str(element.compressionAlgorithm),"None" )

    def testMass(self):
        """ test the mass compression, a positive and negative example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, True, False, 5.*GeV )
        for topo in topos:
            if str(topo)!="[1,1,0][1,0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[e-],[nu]],[[ta+]]]":
                    continue
                masses=element.mother.getMasses()
                dm=(masses[1][1]-masses[1][2])/GeV
                self.assertEqual ( str(element.mother),"[[[e-],[nu]],[[ta+],[ta-]]]")
                self.assertEqual ( str(element.compressionAlgorithm),"mass" )
                self.assertTrue ( dm < 5.0 )

if __name__ == "__main__":
    unittest.main()
