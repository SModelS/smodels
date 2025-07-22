#!/usr/bin/env python3

"""
.. module:: testEvaluationType
   :synopsis: Tests if the evalutionType works as expected (pun intended)

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest

from smodels.statistics.basicStats import NllEvalType, observed, aposteriori, apriori

class EvalTypeTest(unittest.TestCase):

    def testEvalType(self):
        """ test that NllEvalType thing """
        obs = NllEvalType.init ( "observed" )
        obs2 = NllEvalType.init ( False )
        prior = NllEvalType.init ( "prior" )
        prior2 = NllEvalType.init ( True )
        posteriori  =  NllEvalType.init ( "posteriori" )
        posteriori2 =  NllEvalType.init ( "aposteriori" )
        self.assertTrue( obs == obs2 )
        self.assertTrue( obs == obs )
        self.assertTrue( obs == observed )
        self.assertTrue( obs != prior )
        self.assertTrue( obs != posteriori )
        self.assertTrue( prior != posteriori )
        self.assertTrue( prior == prior2 )
        self.assertTrue( prior == apriori )
        self.assertTrue( posteriori == posteriori2 )
        self.assertTrue( posteriori == aposteriori )

if __name__ == "__main__":
    unittest.main()
