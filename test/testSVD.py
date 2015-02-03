#!/usr/bin/env python

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA.
              Will replace the upper limit test

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import math
from smodels.experiment.sVDTrafo import SVDTrafo
from smodels.tools.physicsUnits import GeV, fb

class SVDTest(unittest.TestCase):
    def testSVD(self):
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  3.*fb ], 
                 [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ],  5.*fb ], 
                 [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], 
                 [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], 13.*fb ], 
                 [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], 15.*fb ], 
                 [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], 20.*fb ], 
                 [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ], 8.*fb ], 
                 [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], 10.*fb ], 
                 [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], 12.*fb ], 
                 [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], 15.*fb ], 
                 [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], 17.*fb ], 
                 [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], 19.*fb ], 
                 ]
        trafo=SVDTrafo ( data )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ] ).asNumber(fb), 10. )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ] ).asNumber(fb), 17. )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ] ).asNumber(fb), 11.5 )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,120.*GeV] ] ).asNumber(fb), 11.2 )
        self.assertTrue ( math.isnan ( trafo.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,125.*GeV] ] ).asNumber(fb) ) )
        self.assertTrue ( math.isnan ( trafo.getInterpolatedValue ( [[ 500.*GeV,200.*GeV], [ 500.*GeV,200.*GeV] ] ).asNumber(fb) ) )
        loose=SVDTrafo ( data, accept_errors_upto = 0.05 )
        self.assertAlmostEquals ( loose.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,125.*GeV] ] ).asNumber(fb), 11.35 )
        self.assertTrue ( math.isnan ( loose.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,140.*GeV] ] ).asNumber(fb) ) )

    def testWithEfficiencies(self):
        data = [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .1 ], 
                 [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], .13 ], 
                 [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], .15 ], 
                 [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], .2 ], 
                 [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ], .08 ], 
                 [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], .1 ], 
                 [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], .12 ], 
                 [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], .15 ], 
                 [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], .17 ], 
                 [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], .19 ], 
                 ]
        trafo=SVDTrafo ( data )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ] ), .1 )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ] ), .17 )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,125.*GeV], [ 300.*GeV,125.*GeV] ] ), .115 )
        self.assertAlmostEquals ( trafo.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,120.*GeV] ] ), .112 )
        self.assertTrue ( math.isnan ( trafo.getInterpolatedValue ( [[ 300.*GeV,120.*GeV], [ 300.*GeV,125.*GeV] ] ) ) )
        self.assertTrue ( math.isnan ( trafo.getInterpolatedValue ( [[ 500.*GeV,200.*GeV], [ 500.*GeV,200.*GeV] ] ) ) )

if __name__ == "__main__":
    unittest.main()
