#!/usr/bin/env python

"""
.. module:: regression.testSLHAUniqueName
   :synopsis: Compares the output of the unique naming facility for SLHA files.
    
.. moduleauthor:: Wolfgang magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import setPath
from Tools.VariousHelpers import logging
from Theory import SLHATools


class Test(unittest.TestCase):


    def testUniqueNamingFacility(self):
        slhaFile = "../slha/andrePT1.slha"
        uniqueName = SLHATools.uniqueName ( slhaFile )        
            
        with open('5.log', 'r') as logFile:
            data = logFile.read().replace('\n', '')
        
        self.assertEqual(uniqueName , data)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()