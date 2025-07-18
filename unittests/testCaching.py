#!/usr/bin/env python3

"""
.. module:: testCaching
   :synopsis: Tests the interpolation caching,
        the cache should be slahsed into half after n_stored entries

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.caching import roundCache
from smodels.base.physicsUnits import GeV, fb
from databaseLoader import database
import numpy as np

class CachingTest(unittest.TestCase):

    def testCache(self):

        database.selectExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], 
                    datasetIDs=[None], txnames=["T2bb" ] )
        expRes = database.expResultList
        txname=expRes[0].datasets[0].txnameList[0] # T2bb
        massesvec = [] 
        for i in range(170,290,10):
            massesvec.append ( [ i, 100 ] )
        for masses in massesvec:
            txname.txnameData.getValueFor( np.array([masses,masses]).flatten())

        m = [ 270.0, 100.0, 270.0, 100.0 ]
        result=txname.txnameData.getValueFor(m)
        result = result*txname.y_unit
        self.assertAlmostEqual(result.asNumber(fb) , 459.658) 

    def testRoundCache(self):

        @roundCache(argname='x',argpos=0,digits=1)
        def dummyFunc(x, y):
            return x*y
                        
        self.assertEqual(dummyFunc(1.234,1.0),1.2,3)
        self.assertEqual(dummyFunc(1.234,y=1.0),1.2,3)
        self.assertEqual(dummyFunc(x=1.234,y=1.0),1.2,3)
        self.assertEqual(dummyFunc(1.29,1.0),1.3,3)

        @roundCache(argname='y',argpos=1,digits=2)
        def dummyFunc2(x, y):
            return x*y
        
        self.assertEqual(dummyFunc2(2.,1.2344),2.46,3)
        self.assertEqual(dummyFunc2(x=2.0,y=1.299),2.6,3)


if __name__ == "__main__":
    unittest.main()
