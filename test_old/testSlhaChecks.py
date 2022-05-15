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
        self.assertEqual ( st.status, (1, 'Input file ok') )
        
if __name__ == "__main__":
    unittest.main()
