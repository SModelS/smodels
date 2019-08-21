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
from smodels.tools.caching import Cache
from smodels.tools.physicsUnits import GeV, fb
from databaseLoader import database

class CachingTest(unittest.TestCase):
    def testCache(self):
        Cache.n_stored = 10
        Cache.reset()
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], 
                    datasetIDs=[None], txnames=["T2bb" ] )
        txname=expRes[0].datasets[0].txnameList[0] # T2bb
        massesvec = [] 
        for i in range(170,290,10):
            massesvec.append ( [ i*GeV, 100*GeV ] )
        for masses in massesvec:
            txname.txnameData.getValueFor( [ masses, masses ])
        #    print masses,result,Cache.size()
        self.assertEqual ( Cache.size(), 7 )
        m = [ [ 270*GeV, 100*GeV], [ 270*GeV, 100*GeV ] ]
        result=txname.txnameData.getValueFor(m)
        self.assertAlmostEqual(result.asNumber(fb) , 459.658) 

if __name__ == "__main__":
    unittest.main()
