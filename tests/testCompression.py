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
                # print element,"mother:",element.motherElements,element.compressionAlgorithms,element.weight
                self.assertEqual ( str(element.motherElements[0][1]),"[[[nu],[nu]],[[nu],[ta+]]]" )
                ## 6 neutrino flavors get added!
                self.assertEqual ( len(element.motherElements), 6 ) 
                self.assertEqual ( str(element.motherElements[0][0]),"invisible" )

    def testInvisibleNegative(self):
        """ test the invisible compression, a negative example """
        slhafile="../inputFiles/slha/andrePT4.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[1,1,0][1,1,0]":
                continue
            for e in topo.elementList:
                if str(e)!="[[[nu],[mu-]],[[mu+],[mu-]]]":
                    continue
                #if len(e.motherElements)==1 and e.motherElements[0]=="uncompressed":
                #    print topo,e,e.motherElements
                #self.assertEqual ( str(e.motherElements[0]),"uncompressed" )
                self.assertEqual ( len(e.motherElements),0 )
                #self.assertEqual ( str(e.compressionAlgorithms[0]),"none" )

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
                masses=e.motherElements[0][1].getMasses()
                dm=(masses[1][1]-masses[1][2])/GeV
                self.assertEqual(str(e.motherElements[0][1]),"[[[e-],[nu]],[[ta+],[ta-]]]")
                self.assertEqual(len(e.motherElements),1 )
                self.assertEqual(str(e.motherElements[0][0]),"mass" )
                self.assertTrue ( dm < 5.0 )

if __name__ == "__main__":
    unittest.main()
