#!/usr/bin/env python

"""
.. module:: testXSecComputer
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import set_path
from theory import XSecComputer


class Test(unittest.TestCase):


    def testXSecComputer(self):
        weights = "w= " + XSecComputer.compute(1000, "../slha/andrePT1.slha").weights().__str__()
        with open('1.log', 'r') as logFile:
            data = logFile.read().replace('\n', '')
        self.assertEqual(weights, data)


if __name__ == "__main__":
    unittest.main()
