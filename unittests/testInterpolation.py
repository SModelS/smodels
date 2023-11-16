#!/usr/bin/env python3

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA
              and the triangulation.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.txnameObj import TxNameData
from smodels.base.physicsUnits import GeV, pb, fb
from databaseLoader import database
from unitTestHelpers import theorySMSFromString as fromString
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
import numpy as np

slhafile = './testFiles/slha/lightEWinos.slha'
model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])


class InterpolationTest(unittest.TestCase):
    def testExpected(self):
        expRes = database.getExpResults(analysisIDs=["CMS-PAS-SUS-12-026"], 
                    datasetIDs=[None], txnames=[ "T1tttt" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        sms = fromString('(PV > gluino(1),gluino(2)), (gluino(1) > t+,t-,N1), (gluino(2) > t+,t-,N1)',
                         model=model)
        gluino = model.getParticle(label='gluino')
        n1 = model.getParticle(label='N1')
        gluino.mass = 650.*GeV
        n1.mass = 50.*GeV
        smsMatch = txname.hasSMSas(sms)
        observed = txname.getULFor( smsMatch, expected = False )
        expected = txname.getULFor( smsMatch, expected = True )
        self.assertAlmostEqual(observed.asNumber(fb),49.9,1)
        self.assertAlmostEqual(expected.asNumber(fb),78.5,1)

    def testExpectedFails(self):
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"],
                    datasetIDs=[None], txnames=["T2bb" ] )
        txname=expRes[0].datasets[0].txnameList[0]
        sms = fromString('(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,N1), (sb_1(2) > b,N1)',
                         model=model)
        b1 = model.getParticle(label='sb_1')
        n1 = model.getParticle(label='N1')
        b1.mass = 650.*GeV
        n1.mass = 50.*GeV
        smsMatch = txname.hasSMSas(sms)
        expected = txname.getULFor(smsMatch, expected = True)
        self.assertTrue(expected is None)

    def testInterpolation(self):
        # print database
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"],
                    datasetIDs=[None], txnames=["T2bb" ] )
        #expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T2bb
        result=txname.txnameData.getValueFor([300.,100.,300.,100.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb),0.162457, places=4 )
        result=txname.txnameData.getValueFor([300.,125.,300.,125.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb),0.237745, places=4 )

    def test6D(self):
        # print database
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"],
                txnames=[ "T6bbWW" ] )
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T6bbWW
        result=txname.txnameData.getValueFor([ 300.,105.,100.,300.,105.,100.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb),0.176266, places=4 )
        result=txname.txnameData.getValueFor([ 300.,270.,200.,300.,270.,200.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb), 87.0403, places=3 )
        result=txname.txnameData.getValueFor([300.,270.,200.,300.,271.,200.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb), 88.6505675, places=3 )


    def testOutsidePlane(self):
        expRes = database.getExpResults( analysisIDs=["ATLAS-SUSY-2013-05"],
                                         txnames=["T2bb" ] )
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname=expRes[0].datasets[0].txnameList[0] # T6bbWW
        result=txname.txnameData.getValueFor([300.,127.,300.,127.5])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb),0.24452092,places=4 )
        result=txname.txnameData.getValueFor([600.,120.,600.,130.])
        result = result*txname.y_unit
        self.assertAlmostEqual( result.asNumber(pb),0.0197154,places=4 )
        result=txname.txnameData.getValueFor([300.,120.,300.,130.])
        self.assertTrue ( result is None )

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

        xvalues = []
        yvalues = []
        for pt in data:
            mass = np.array(pt[0]).flatten()
            xvalues.append([m.asNumber(GeV) for m in mass])
            yvalues.append(pt[1].asNumber(pb))
        txnameData=TxNameData (x=xvalues,y=yvalues, txdataId='test')
        result=txnameData.getValueFor([300.,125.,300.,125.])
        self.assertAlmostEqual(result,0.0115, 3)

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

        xvalues = []
        yvalues = []
        for pt in data:
            mass = np.array(pt[0]).flatten()
            xvalues.append([m.asNumber(GeV) for m in mass])
            yvalues.append(pt[1])
        txnameData=TxNameData (x=xvalues,y=yvalues, txdataId='test2')
        result=txnameData.getValueFor([300.,125.,300.,123.])
        self.assertAlmostEqual(result,0.1144 )

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
        xvalues = []
        yvalues = []
        for pt in data:
            mass = np.array(pt[0]).flatten()
            xvalues.append([m.asNumber(GeV) for m in mass])
            yvalues.append(pt[1])
        txnameData=TxNameData (x=xvalues,y=yvalues, txdataId='test3')
        result=txnameData.getValueFor([300.,125.,300.,100.])
        self.assertEqual(result, None)
       

if __name__ == "__main__":
    unittest.main()
