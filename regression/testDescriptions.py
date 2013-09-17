#!/usr/bin/env python

"""
.. module:: testDescriptions
   :synopsis: Compare topology descriptions with given values.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import setPath
from Experiment import TxNames
from Theory import LHEReader, TopologyBuilder


class Test(unittest.TestCase):


    def testDescriptions(self):
        topologyList = ['T1','T2','T1tttt', 'T2tt','T3W', 'T5WW', 'TChiWZ', 'T1bbbb', 'T2bb', 'T5WZ', 'T3Wb', 'T3Z', 'T5ZZ', 'T6bbZZ', 'TChiWW', 'TSlepSlep']
        for topology in topologyList:
            event = LHEReader.LHEReader("../lhe/%s_1.lhe" % topology).next()
            smsTopology = TopologyBuilder.fromEvent(event, {})
            result = TxNames.getTx(smsTopology[0].leadingElement())
            self.assertEqual(result, topology)


if __name__ == "__main__":
    unittest.main()