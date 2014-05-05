#!/usr/bin/env python

"""
.. module:: testXSecComputer
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest
import setPath
from smodels.tools import xsecComputer
import logging.config

class Test(unittest.TestCase):
    def TestXSecComputer(self):
        w = xsecComputer.computeXSec(8,0,1000, "../inputFiles/slha/andrePT4.slha").getDictionary()
        w8lo= 1000 * w[(1000023, 1000024)]['8 TeV (LO)'].asNumber() 
        self.assertAlmostEqual(w8lo, 35.014621117 )  ## 35.01 fb

    def testSecond (self):
        logging.config.fileConfig('./logging.conf' )
        w = xsecComputer.computeXSec(8,2,1000, "../inputFiles/slha/andrePT4.slha").getDictionary()
        print "w=",w
        self.assertEqual(3*4,12)

if __name__ == "__main__":
    unittest.main()
