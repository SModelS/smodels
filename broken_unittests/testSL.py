#!/usr/bin/env python3

"""
.. module:: testSL
   :synopsis: Test the Simplified Likelihoods

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys

sys.path.insert(0, "../")
import unittest
from smodels.tools.simplifiedLikelihoods import Data, UpperLimitComputer, \
      LikelihoodComputer
from smodels.tools.physicsUnits import fb
from numpy import sqrt
import numpy as np


class SLTest(unittest.TestCase):
    def testPathologicalModel(self):
        C = [1.0]
        m = Data(
            observed=[0],
            backgrounds=[0.0],
            covariance=C,
            third_moment=[0.0] * 8,
            nsignal=[x / 100.0 for x in [0.0]],
            name="pathological model",
            deltas_rel=0.0,
        )
        m.zeroSignal()

        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ul = ulComp.getUpperLimitOnMu(m, marginalize=True )
        ulProf = ulComp.getUpperLimitOnMu ( m, marginalize=False )
        self.assertEqual(ul, None )
        self.assertEqual(ulProf, None )

    def testPathologicalModel2(self):
        C = [1.0]
        m = Data(
            observed=[0],
            backgrounds=[0.0],
            covariance=C,
            third_moment=[0.0] * 8,
            nsignal=[x / 100.0 for x in [0.1]],
            name="pathological model 2",
            deltas_rel=0.0,
        )
        m.zeroSignal()

        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ul = ulComp.getUpperLimitOnMu(m, marginalize=True )
        ulProf = ulComp.getUpperLimitOnMu(m, marginalize=False )
        self.assertAlmostEqual(ul, 3049.617040116002, 1 )
        self.assertAlmostEqual(ulProf, 1920.7547785412298, 1 )

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
        ulProf = ulComp.getUpperLimitOnMu ( m, marginalize=False )
        self.assertAlmostEqual( ulProf / 131.58637474312224, 1.0, 2 )
        ul = ulComp.getUpperLimitOnMu ( m, marginalize = True )
        self.assertAlmostEqual( ul / 132.75018479789006, 1., 1 )

    def createModel(self,n=3):
        import model_90 as m9

        S = m9.third_moment.tolist()[:n]
        D = m9.observed.tolist()[:n]
        B = m9.background.tolist()[:n]
        sig = [x / 100.0 for x in m9.signal.tolist()[:n]]
        C_ = m9.covariance.tolist()
        ncov = int(sqrt(len(C_)))
        C = []
        for i in range(n):
            C.append(C_[ncov * i : ncov * i + n])
        m = Data( observed=D, backgrounds=B, covariance=C, third_moment=S,
            nsignal=sig, name="model%d" % n, deltas_rel=0.0,)
        return m

    def testModel3(self):

        """ take first n SRs of model-90 """
        m = self.createModel ( 3 )
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        lComp = LikelihoodComputer( m )
        ulProf = ulComp.getUpperLimitOnMu( m, marginalize=False )
        self.assertAlmostEqual( ulProf / 2168.8056715301045, 1.0, 3 )
        ul = ulComp.getUpperLimitOnMu( m, marginalize=True )
        ## Nick's profiling code gets for n=3 ul=2135.66
        self.assertAlmostEqual(ul / 2195.3529619253704, 1.0, 1)
        lmax = lComp.lmax ( )
        self.assertAlmostEqual( lmax, 2.1722054152e-09, 7 )
        self.assertAlmostEqual( lComp.muhat, 0. )
        # self.assertAlmostEqual( lComp.sigma_mu, 798.9887891147, 2 )
        self.assertAlmostEqual( lComp.sigma_mu, 800.1380826235132, 2 )
        lmax = lComp.lmax ( allowNegativeSignals=True )
        self.assertAlmostEqual( lmax, 2.1708552631256182e-09, 12 )
        #self.assertAlmostEqual( lComp.muhat, -71.523083468, 7 )
        self.assertAlmostEqual( lComp.muhat, -72.63852360245156, 7 )
        #self.assertAlmostEqual( lComp.sigma_mu, 795.0121298843319 )
        self.assertAlmostEqual( lComp.sigma_mu, 796.122901385052, 4 )

    def xestModel1(self):
        """ take first SR of model-90 """
        m = self.createModel ( 1 )
        m.nsignal[0]=1.
        lComp = LikelihoodComputer( m )
        lmax = lComp.lmax ( m.nsignal )
        self.assertAlmostEqual( lmax, 0.0003441355122238784 )
        self.assertAlmostEqual( lComp.muhat, 1. )
        self.assertAlmostEqual( lComp.sigma_mu, 32.31764780503341 )
        lmax = lComp.lmax ( m.nsignal, allowNegativeSignals=True )
        self.assertAlmostEqual( lmax, 0.0003441355122238784 )
        self.assertAlmostEqual( lComp.muhat, 1. )
        self.assertAlmostEqual( lComp.sigma_mu, 32.31764780503341 )

        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        #ulProf = ulComp.getUpperLimitOnMu( m, marginalize=False )
        #self.assertAlmostEqual( ulProf/54.793636190198924, 1.0, 3 )
        ul = ulComp.getUpperLimitOnMu( m, marginalize=False )
        ## Nick's profiling code gets for n=3 ul=2135.66
        self.assertAlmostEqual(ul / 61.26914, 1.0, 1)

    def testModel10(self):

        """ take first 10 SRs of model-90 """
        m = self.createModel ( 10 )
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ulProf = ulComp.getUpperLimitOnMu( m, marginalize=False )
        self.assertAlmostEqual( ulProf / 365.6091713369213, 1.0, 2 )
        ul = ulComp.getUpperLimitOnMu(m,marginalize=True)
        ## Nick's profiling code gets for n=10 ul=357.568
        self.assertAlmostEqual(ul / 371.047747734418, 1.0, 1)

    def testModel40(self):
        m = self.createModel(40)
        import time

        ulComp = UpperLimitComputer(ntoys=20000, cl=.95 )
        ulProf = ulComp.getUpperLimitOnMu ( m, marginalize=False )
        self.assertAlmostEqual( ulProf / 61.53473539725907, 1.0, 2 )
        ul = ulComp.getUpperLimitOnMu ( m, marginalize=True )
        self.assertAlmostEqual ( ul/62.943195310667136, 1., 1 )

    def testTrivialModel ( self ):
        def pprint ( *args ):
            return
            print ( " ".join ( map ( str, args ) ) )
        name = "Trivial Model"
        observed = [ 3., 3. ]
        background = [ 1., 1. ]
        covariance = [[ 1., 0. ], [ 0., 1. ] ]
        signal = [ 2., 2. ]
        pprint ( "nsignal", sum(signal) )
        m=Data ( observed=observed,
                  backgrounds=background,
                  covariance= covariance,
                  third_moment = [ 0. ] * 2,
                  nsignal= signal,
                  name="trivial model",deltas_rel=0., lumi = 100./fb )
        lComp = LikelihoodComputer( m )
        theta_hat = lComp.findThetaHat( signal )
        pprint ( "theta_hat", theta_hat )
        lmax = lComp.lmax()
        lm = lComp.likelihood ( 1. )
        muhat = lComp.muhat
        sigma_mu = lComp.sigma_mu
        pprint ( "muhat", muhat, "sigma_mu", sigma_mu )
        theta_hat = lComp.findThetaHat( muhat * np.array ( signal ) )
        pprint ( "theta_hat", theta_hat )
        # mu_hat is (observed - background) / signal
        # (as it's the same for all regions)
        # it's 1.0
        ulComp = UpperLimitComputer(ntoys=10000, cl=.95 )
        ulProf = ulComp.getUpperLimitOnMu ( m, marginalize=False )
        pprint ( "ulProf", ulProf )
        self.assertAlmostEqual ( lmax, lm, 3 )
        self.assertAlmostEqual ( lComp.muhat, 1., 5 )
        self.assertAlmostEqual ( lComp.sigma_mu, 0.655333438756686, 5 )
        # self.assertAlmostEqual ( lComp.sigma_mu, 0.6123724356957945 )
        self.assertAlmostEqual ( ulProf, 2.4412989438119546, 5 )
        # self.assertAlmostEqual ( ulProf, 2.537934980801342, 5 )


if __name__ == "__main__":
    unittest.main()
