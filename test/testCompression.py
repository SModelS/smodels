#!/usr/bin/env python

"""
.. module:: testCompression
   :synopsis: Checks the compression algorithms
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import GeV, fb
import unittest

class CompressionTest(unittest.TestCase):
    from smodels.tools.smodelsLogging import logger

    def testInvisiblePositive(self):
        """ test the invisible compression, a positive example """
        slhafile="../inputFiles/slha/higgsinoStop.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[0][0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[],[]]":
                    continue
                #print "m00=",str(element.motherElements[0][0])
                self.assertEqual ( str(element.motherElements[0][1]),"[[],[[nu,nu]]]")
                ## 3 neutrino flavors get added!
                self.assertEqual ( len(element.motherElements), 3 ) 
                self.assertEqual ( str(element.motherElements[0][0]),"invisible" )

    def testInvisibleNegative(self):
        """ test the invisible compression, a negative example """
        slhafile="../inputFiles/slha/higgsinoStop.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, True, 5.*GeV )
        for topo in topos:
            if str(topo)!="[1,1,0][1,1,0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[jet],[W+]],[[t-],[t+]]]":
                    continue
                #print
                #print topo,element,"mother:",len(element.motherElements),element.motherElements
                #for x in element.motherElements: 
                #    print "m0",str(x[0]),str(x[1])
                #if len(e.motherElements)==1 and e.motherElements[0]=="uncompressed":
                #    print topo,e,e.motherElements
                #self.assertEqual ( str(e.motherElements[0]),"uncompressed" )
                self.assertEqual ( len(element.motherElements),0 )
                #self.assertEqual ( str(element.motherElements[0][1]),"[]")
                #self.assertEqual ( str(e.compressionAlgorithms[0]),"none" )

    def testMass(self):
        """ test the mass compression, a positive example """
        slhafile="../inputFiles/slha/higgsinoStop.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, True, False, 5.*GeV )
        for topo in topos:
            if str(topo)!="[1,0][1,0]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[b]],[[b]]]":
                    continue
                masses=element.motherElements[0][1].getMasses()
                #for x in element.motherElements:
                #    if str(x[0])=="mass":
                #        print str(topo),str(element),str(x[1])
                #        print "    --- ",masses
                dm=abs(masses[0][1]-masses[0][2])/GeV
                ## #self.assertEqual(str(element.motherElements[0][1]),"[[[e-],[nu]],[[ta+],[ta-]]]")
                self.assertEqual(len(element.motherElements),25 )
                self.assertEqual(str(element.motherElements[0][0]),"mass" )
                self.assertTrue ( dm < 5.0 )

if __name__ == "__main__":
    unittest.main()
