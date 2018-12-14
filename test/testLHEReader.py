#!/usr/bin/env python3

"""
.. module:: testAsciiGraph
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import unittest
import sys
sys.path.insert(0,"../")
from smodels.theory import lheReader
import pyslha
import logging as logger


def compareDecays(decaysA,decaysB,allowDiff=0.01,minBR=1e-5):
    """
    Compares two lists of pyslha.Decay objects. Each list should correpond
    to the same parent particle.
    """
    
    dA = [decay for decay in decaysA if decay.br > minBR]
    dB = [decay for decay in decaysB if decay.br > minBR]
    #First sort both decays by final state
    dA = sorted(dA, key=lambda x: sorted(x.ids))
    dB = sorted(dB, key=lambda x: sorted(x.ids))        

    if len(dA) != len(dB):
        logger.error("Decay length differ")
        
        print(dA)
        print(dB)
        return False


    for idec,decayA in enumerate(dA):
        decayB = dB[idec]
        if sorted(decayA.ids) != sorted(decayB.ids):
            logger.error("Decay final states differ: %s and %s" %(str(decayA),str(decayB)))
            return False
        
        if 2*abs(decayA.br-decayB.br)/abs(decayA.br+decayB.br) > allowDiff:
            logger.error("Decay BRs differ: %s and %s" %(str(decayA),str(decayB)))
            return False
        
    return True

class LHEReaderTest(unittest.TestCase):

    def testReader(self):
        
        filename = "./testFiles/lhe/simplyGluino.lhe"
        reader = lheReader.LheReader(filename)
        events = [event for event in reader]
        reader.close()
        self.assertEqual(len(events),5)
        

    def testEventDictionaries(self):

        filename = "./testFiles/lhe/simplyGluino.lhe"
        reader = lheReader.LheReader(filename)
        events = [event for event in reader]
        reader.close()
        
        massDict,decayDict = lheReader.getDictionariesFromEvent(events[0])
        self.assertEqual(massDict,{-1: [0.33], 1000021: [675.0], 1000022: [200.0], 1: [0.33]})
        gluinoDec = [pyslha.Decay(br=1.,nda=3,ids=[-1,1,1000022],parentid=1000021)]*2
        self.assertTrue(compareDecays(gluinoDec,decayDict[1000021]))
        self.assertEqual(decayDict[1000022],[])
        
    def testDictionaries(self):

        filename = "./testFiles/lhe/simplyGluino.lhe"
        massDict,decayDict = lheReader.getDictionariesFrom(filename)
        self.assertEqual(massDict,{1000021: 675.0, 1000022: 200.0, 1: 0.33, 2 : 0.33})
        gluinoDecs = [pyslha.Decay(br=0.3,nda=3,ids=[-1,1,1000022],parentid=1000021),
                      pyslha.Decay(br=0.7,nda=3,ids=[-2,2,1000022],parentid=1000021)]
        self.assertEqual(len(decayDict[1000021].decays),len(gluinoDecs))
        self.assertTrue(compareDecays(gluinoDecs,decayDict[1000021].decays))            
        
        re = pyslha.readSLHAFile("./testFiles/slha/gluino_squarks.slha")

        filename = "./testFiles/lhe/gluino_squarks.lhe"
        massDict,decayDict = lheReader.getDictionariesFrom(filename)
        for pdg in massDict:
            if pdg < 100000:
                continue            
            self.assertAlmostEqual(re.blocks['MASS'][pdg],massDict[pdg])
            
        #Expected answer:
        decayRes = {}
        decayRes[1000024] = [pyslha.Decay(br=1.0000,ids=[1000022, 24],parentid=1000024,nda=2)]
        decayRes[1000023] = [pyslha.Decay(br=0.8571,ids=[1000022, 25],parentid=1000023,nda=2),
                             pyslha.Decay(br=0.1429,ids=[1000022, 23],parentid=1000023,nda=2)]
        decayRes[1000001] = [pyslha.Decay(br=1.0000,ids=[1000021, 1],parentid=1000001,nda=2)]
        decayRes[1000002] = [pyslha.Decay(br=0.5000,ids=[1000024, 1],parentid=1000002,nda=2),
                             pyslha.Decay(br=0.2500,ids=[1000021, 2],parentid=1000002,nda=2),
                             pyslha.Decay(br=0.2500,ids=[1000023, 2],parentid=1000002,nda=2)]
        decayRes[2000002] = [pyslha.Decay(br=1.0000,ids=[1000021, 2],parentid=2000002,nda=2)]
        decayRes[1000021] = [pyslha.Decay(br=0.2500,ids=[-1000024, 4, -3],parentid=1000021,nda=3),
                             pyslha.Decay(br=0.3750,ids=[1000024, -2, 1],parentid=1000021,nda=3),
                             pyslha.Decay(br=0.1250,ids=[1000024, -4, 3],parentid=1000021,nda=3),
                             pyslha.Decay(br=0.1250,ids=[1000023, -2, 2],parentid=1000021,nda=3),
                             pyslha.Decay(br=0.1250,ids=[1000023, -6, 6],parentid=1000021,nda=3)]

        for pdg in decayRes:
            self.assertTrue(compareDecays(decayDict[pdg].decays,decayRes[pdg]))
            self.assertEqual(decayDict[pdg].totalwidth,float('inf'))
            
        for pdg in decayDict:
            if not decayDict[pdg].decays:
                self.assertEqual(decayDict[pdg].totalwidth,0.)

if __name__ == "__main__":
    unittest.main()
