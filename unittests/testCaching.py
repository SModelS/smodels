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

        # self.assertEqual ( txname.txnameData.getValueFor.cache_parameters(), {'maxsize': 128, 'typed': False} )
        m = [ 270.0, 100.0, 270.0, 100.0 ]
        for i in range(3):
            result=txname.txnameData.getValueFor(m)
        result = result*txname.y_unit
        self.assertAlmostEqual(result.asNumber(fb) , 459.658) 

        self.assertEqual ( txname.txnameData.getValueFor.__wrapped__.cache_info().hits, 3 )

if __name__ == "__main__":
    unittest.main()
