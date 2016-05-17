#!/usr/bin/env python

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA.
              Will replace the upper limit test

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.databaseObj import Database
from smodels.experiment.txnameObj import TxNameData
from smodels.tools.physicsUnits import GeV, TeV, pb, fb
import math

class InterpolationTest(unittest.TestCase):
    def testInterpolation(self):
        self.database = Database ( "./database/database.pcl" )
        # print database
        expRes = self.database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], 
                    datasetIDs=[None], txnames=["T2bb" ] )
        #expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T2bb
        result=txname.txnameData.getValueFor([[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.162457 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.237745 )
    def test6D(self):
        self.database = Database ( "./database/database.pcl" )
        # print database
        expRes = self.database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=[ "T6bbWW" ] )
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T6bbWW
        result=txname.txnameData.getValueFor([[ 300.*GeV,105.*GeV,100.*GeV], [ 300.*GeV,105.*GeV,100.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.176266 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,270.*GeV,200.*GeV], [ 300.*GeV,270.*GeV,200.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb), 87.0403 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,270.*GeV,200.*GeV], [ 300.*GeV,271.*GeV,200.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb), 88.6505675 )
    def testOutsidePlane(self):
        self.database = Database ( "./database/database.pcl" )
        expRes = self.database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=["T2bb" ] )
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T6bbWW
        result=txname.txnameData.getValueFor([[ 300.*GeV,127.*GeV], [ 300.*GeV,127.5*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.24452092000000109 )
        result=txname.txnameData.getValueFor([[ 600.*GeV,120.*GeV], [ 600.*GeV,130.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.0197154 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,120.*GeV], [ 300.*GeV,130.*GeV] ])
        self.assertTrue ( result == None )

    def testWithDirectData(self):
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  3.*fb ], 
             [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ],  5.*fb ], 
             [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], 
             [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], 13.*fb ], 
             [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], 15.*fb ], 
             [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], 20.*fb ], 
             [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ],  8.*fb ], 
             [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], 10.*fb ], 
             [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], 12.*fb ], 
             [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], 15.*fb ], 
             [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], 17.*fb ], 
             [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], 19.*fb ], ]
        txnameData=TxNameData ( data ) ## "upperlimit", data )
        result=txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.0115 ) 

    def testEfficiencyMaps(self):
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  .03 ], 
             [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ], .05 ], 
             [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .10 ], 
             [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], .13 ], 
             [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], .15 ], 
             [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], .20 ], 
             [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ], .08 ], 
             [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], .10 ], 
             [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], .12 ], 
             [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], .15 ], 
             [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], .17 ], 
             [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], .19 ], ]
        txnameData=TxNameData ( data) ## "upperlimit", data )
        result=txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ])
        self.assertAlmostEquals( result,0.115 ) 
        
    def testOutsideConvexHull( self ):
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  .03 ], 
             [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ], .05 ], 
             [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .10 ], 
             [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], .13 ], 
             [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], .15 ], 
             [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], .20 ], 
             [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ], .08 ], 
             [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], .10 ], 
             [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], .12 ], 
             [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], .15 ], 
             [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], .17 ], 
             [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], .19 ], ]
        txnameData=TxNameData ( data)
        result=txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,123.*GeV] ])
        self.assertAlmostEquals( result,0.1144 ) 

    def testOutsideConvexHull2( self ):
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  .03 ], 
             [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ], .05 ], 
             [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .10 ], 
             [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], .13 ], 
             [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], .15 ], 
             [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], .20 ], 
             [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ], .08 ], 
             [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], .10 ], 
             [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], .12 ], 
             [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], .15 ], 
             [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], .17 ], 
             [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], .19 ], ]
        txnameData=TxNameData ( data)
        result=txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,100.*GeV] ])
        self.assertEquals ( result, None )
       

if __name__ == "__main__":
    unittest.main()
