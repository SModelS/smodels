#!/usr/bin/env python

"""
.. module:: testCaching
   :synopsis: Tests the interpolation caching,
        the cache should be slahsed into half after n_stored entries

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.installation import installDirectory as idir
from smodels.tools.caching import Cache
from smodels.tools.physicsUnits import GeV, fb
from databaseLoader import database
from os.path import join

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
            result=txname.txnameData.getValueFor( [ masses, masses ])
        #    print masses,result,Cache.size()
        self.assertEquals ( Cache.size(), 7 )
        m = [ [ 270*GeV, 100*GeV], [ 270*GeV, 100*GeV ] ]
        result=txname.txnameData.getValueFor( m )
        self.assertAlmostEquals ( result.asNumber ( fb ) , 459.658 ) 

if __name__ == "__main__":
    unittest.main()
