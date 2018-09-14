#!/usr/bin/env python3

"""
.. module:: testSlhaChecks
   :synopsis: Tests the slha checker

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools import slhaChecks

class SlhaTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "./testFiles/slha/gluino_squarks.slha"
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals ( st.status, (1, 'Input file ok') )
        
    def testBadFile(self):
        filename = "./testFiles/slha/nobdecay.slha"
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals (st.status, (-1, '#ERROR: special signatures in this point.\n#Warnings:\n##Visible decays of longlived particles / stable charged particles: [1000005]\n#1000005 : c*tau = inf\n\n'))

if __name__ == "__main__":
    unittest.main()
