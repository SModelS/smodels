#!/usr/bin/env python

"""
.. module:: testAsciiGraph
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.installation import installDirectory
from smodels.theory import lheReader
import pyslha

class LHEReaderTest(unittest.TestCase):

    def testReader(self):
        
        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory())
        reader = lheReader.LheReader(filename)
        events = [event for event in reader]
        self.assertEqual(len(events),5)
        

    def testEventDictionaries(self):

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory())
        reader = lheReader.LheReader(filename)
        events = [event for event in reader]
        
        
        massDict,decayDict = lheReader.getDictionariesFromEvent(events[0])
        self.assertEqual(massDict,{-1: [0.33], 1000021: [675.0], 1000022: [200.0], 1: [0.33]})
        gluinoDec = pyslha.Decay(br=1.,nda=3,ids=[-1,1,1000022],parentid=1000021)
        for decay in decayDict[1000021]:
            self.assertEqual(sorted(decay.ids), sorted(gluinoDec.ids))
            self.assertEqual(decay.parentid, gluinoDec.parentid)
            self.assertEqual(decay.br, gluinoDec.br)

        self.assertEqual(decayDict[1000022],[])
        
    def testDictionaries(self):

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory())
        massDict,decayDict = lheReader.getDictionariesFrom(filename)
        self.assertEqual(massDict,{1000021: 675.0, 1000022: 200.0, 1: 0.33, 2 : 0.33})
        gluinoDecs = [pyslha.Decay(br=0.3,nda=3,ids=[-1,1,1000022],parentid=1000021),
                      pyslha.Decay(br=0.7,nda=3,ids=[-2,2,1000022],parentid=1000021)]
        self.assertEqual(len(decayDict[1000021]),len(gluinoDecs))
        for dec in decayDict[1000021]:
            for decB in gluinoDecs:
                if sorted(dec.ids) == sorted(decB.ids):
                    self.assertAlmostEqual(dec.br,decB.br)
                    self.assertAlmostEqual(dec.nda,decB.nda)
                    

if __name__ == "__main__":
    unittest.main()
