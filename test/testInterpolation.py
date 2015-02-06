#!/usr/bin/env python

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA.
              Will replace the upper limit test

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
from smodels.experiment.databaseObjects import DataBase
from smodels.tools.physicsUnits import GeV, TeV, pb
import math

class InterpolationTest(unittest.TestCase):
    def testInterpolation(self):
        self.database = DataBase ( "./database/" )
        # print database
        listOfExpRes = self.database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=["T2bb","T6bbWW" ] )
        expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes.txnames[1] # T2bb
        result=txname.txnameData.getValueFor([[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.162457 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.237745 )
    def testOutsidePlane(self):
        self.database = DataBase ( "./database/" )
        listOfExpRes = self.database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=["T2bb","T6bbWW" ] )
        expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes.txnames[1] # T2bb
        result=txname.txnameData.getValueFor([[ 300.*GeV,126.*GeV], [ 300.*GeV,128.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.24376804 )
        result=txname.txnameData.getValueFor([[ 600.*GeV,120.*GeV], [ 600.*GeV,130.*GeV] ])
        self.assertAlmostEquals( result.asNumber(pb),0.0197154 )
        result=txname.txnameData.getValueFor([[ 300.*GeV,120.*GeV], [ 300.*GeV,130.*GeV] ])
        self.assertTrue ( math.isnan ( result.asNumber(pb) ) )

if __name__ == "__main__":
    unittest.main()
