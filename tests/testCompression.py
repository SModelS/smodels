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
                print element,"mother:",element.mothers[0],element.compressionAlgorithms
                self.assertEqual ( str(element.mothers[0]),"[[[nu],[nu]],[[nu],[ta+]]]" )
                self.assertEqual ( len(element.mothers), 1 )
                self.assertEqual ( str(element.compressionAlgorithms[0]),"invisible" )
                self.assertEqual ( len(element.compressionAlgorithms), 1 )

    def testInvisibleNegative(self):
        """ test the invisible compression, a negative example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[0][1,1,0]":
                continue
            for e in topo.elementList:
                if str(e)!="[[],[[nu],[mu-]]]":
                    continue
                ## print element,"mother:",element.mothers,element.compressionAlgorithm
                self.assertEqual ( str(e.mothers[0]),"self" )
                self.assertEqual ( len(e.mothers),1 )
                self.assertEqual ( str(e.compressionAlgorithms[0]),"direct" )

    def testMass(self):
        """ test the mass compression, a positive and negative example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, True, False, 5.*GeV )
        for topo in topos:
            if str(topo)!="[1,1,0][1,0]":
                continue
            for e in topo.elementList:
                if str(e)!="[[[e-],[nu]],[[ta+]]]":
                    continue
                masses=e.mothers[0].getMasses()
                dm=(masses[1][1]-masses[1][2])/GeV
                self.assertEqual(str(e.mothers[0]),"[[[e-],[nu]],[[ta+],[ta-]]]")
                self.assertEqual(len(e.mothers),1 )
                self.assertEqual(str(e.compressionAlgorithms[0]),"mass" )
                self.assertTrue ( dm < 5.0 )

if __name__ == "__main__":
    unittest.main()
