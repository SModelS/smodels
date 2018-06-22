#!/usr/bin/env python3

"""
.. module:: testSL
   :synopsis: Test the Simplified Likelihoods

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.simplifiedLikelihoods import Data, UpperLimitComputer
from numpy  import sqrt

class SLTest(unittest.TestCase):

    def testPathologicalModel(self):
        C=[ 1. ]
        m=Data ( observed=[0],
                  backgrounds=[.0],
                  covariance= C,
                  third_moment = [ 0. ] * 8,
                  nsignal=[x/100. for x in [0.] ],
                  name="pathological model",deltas_rel=0. )
        m.zeroSignal()
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ul = ulComp.ulSigma(m, marginalize=True )
        ulProf = ulComp.ulSigma ( m, marginalize=False )
        self.assertEqual(ul, None )
        self.assertEqual(ulProf, None )

    def testPathologicalModel2(self):
        C=[ 1. ]
        m=Data(observed=[0],
                  backgrounds=[.0],
                  covariance= C,
                  third_moment = [ 0. ] * 8,
                  nsignal=[x/100. for x in [0.1] ],
                  name="pathological model 2",deltas_rel=0. )
        m.zeroSignal()
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ul = ulComp.ulSigma(m, marginalize=True )
        ulProf = ulComp.ulSigma(m, marginalize=False )
        self.assertAlmostEqual(ul/(3049.*sum(m.nsignal)), 1., 1 )
        self.assertAlmostEqual(ulProf/(1920.*sum(m.nsignal)), 1., 1 )

    def testModel8(self):
        C=[ 18774.2, -2866.97,-5807.3,-4460.52,-2777.25,-1572.97, -846.653, -442.531,
           -2866.97, 496.273, 900.195, 667.591, 403.92, 222.614, 116.779, 59.5958, 
           -5807.3, 900.195, 1799.56, 1376.77, 854.448, 482.435, 258.92, 134.975, 
           -4460.52, 667.591, 1376.77, 1063.03, 664.527, 377.714, 203.967, 106.926, 
           -2777.25, 403.92, 854.448, 664.527, 417.837, 238.76, 129.55, 68.2075, 
           -1572.97, 222.614, 482.435, 377.714, 238.76, 137.151, 74.7665, 39.5247,
           -846.653, 116.779, 258.92, 203.967, 129.55, 74.7665, 40.9423, 21.7285, 
           -442.531, 59.5958, 134.975, 106.926, 68.2075, 39.5247, 21.7285, 11.5732]
        nsignal = [ x/100. for x in [47,29.4,21.1,14.3,9.4,7.1,4.7,4.3] ]
        m=Data ( observed=[1964,877,354,182,82,36,15,11],
                  backgrounds=[2006.4,836.4,350.,147.1,62.0,26.2,11.1,4.7],
                  covariance= C,
                  third_moment = [ 0. ] * 8,
                  nsignal= nsignal,
                  name="CMS-NOTE-2017-001 model",deltas_rel=0. )
        ulComp = UpperLimitComputer (ntoys=2000, cl=.95 )
        ul = ulComp.ulSigma ( m )
        self.assertAlmostEqual( ul/(131.828*sum(nsignal)), 1.0, 1 )
        ulProf = ulComp.ulSigma ( m, marginalize=False )
        #print ( "ul,ulprof=", ul,ulProf )
        self.assertAlmostEqual( ulProf/(131.828*sum(nsignal)), 1.0, 1 )

    def createModel(self,n=3):
        import model_90 as m9
        S=m9.third_moment.tolist()[:n]
        D=m9.observed.tolist()[:n]
        B=m9.background.tolist()[:n]
        sig=[ x/100. for x in m9.signal.tolist()[:n] ]
        C_=m9.covariance.tolist()
        ncov=int(sqrt(len(C_)))
        C=[]
        for i in range(n):
            C.append ( C_[ncov*i:ncov*i+n] )
        m = Data ( observed=D, backgrounds=B, covariance=C, third_moment=S, 
                    nsignal=sig, name="model%d" % n, deltas_rel=0. )
        return m

    def testModel3(self):
        """ take first n SRs of model-90 """
        m = self.createModel ( 3 )
        import time
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        t0=time.time()
        ul = ulComp.ulSigma( m )
        t1=time.time()
        #print ( "ul=%s, t=%s" % ( ul, t1-t0 ) )
        ## Nick's profiling code gets for n=3 ul=2135.66
        self.assertAlmostEqual( ul/(2135.66*sum(m.nsignal)), 1.0, 1 )
        ulProf = ulComp.ulSigma( m, marginalize=False )
        t2=time.time()
        #print ( "ulProf,t=", ulProf, t2-t1 )
        self.assertAlmostEqual( ulProf/(2135.66*sum(m.nsignal)), 1.0, 1 )

    def testModel10(self):
        """ take first 10 SRs of model-90 """
        m = self.createModel ( 10 )
        import time
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        t0=time.time()
        ul = ulComp.ulSigma(m)
        t1=time.time()
        #print ( "ul=%s, t=%s" % ( ul, t1-t0 ) )
        ## Nick's profiling code gets for n=10 ul=357.568
        self.assertAlmostEqual(ul/(357.*sum(m.nsignal)), 1.0, 1 )
        ulProf = ulComp.ulSigma( m, marginalize=False )
        t2=time.time()
        #print ( "ulProf,t=", ulProf, t2-t1 )
        self.assertAlmostEqual( ulProf/(350.*sum(m.nsignal)), 1.0, 1 )

    def testModel40(self):
        m = self.createModel ( 40 )
        import time
        ulComp = UpperLimitComputer(ntoys=5000, cl=.95 )
        ul = ulComp.ulSigma ( m )
        ulProf = ulComp.ulSigma ( m, marginalize=False )
        self.assertAlmostEqual ( ul/(66.*sum(m.nsignal)), 1., 1 )
        self.assertAlmostEqual( ulProf/(63.*sum(m.nsignal)), 1.0, 1 )

if __name__ == "__main__":
    unittest.main()
