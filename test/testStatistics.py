#!/usr/bin/env python

"""
.. module:: testStatistics
   :synopsis: Tests the statistics module.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
from smodels.tools import statistics
import scipy.stats
import math

class StatisticsTest(unittest.TestCase):
#    def testCLInterval(self):
#        re = statistics.computeCLInterval ( 100., 100., 1. )
#        self.assertAlmostEqual ( re, 18.0792727821 )
#
#    def testBayesianLimit(self):
#        re = statistics.bayesianUpperLimit ( 100, 0., 100., 0. )
#        self.assertAlmostEqual ( re, 21.4256127382 )
#
#    def testUL(self):
#        re = statistics.getUL ( 100, 100., 0. )
#
    def testCoverage(self):
        lambdaBG=20
        lambdaSig=5
        nBG=scipy.stats.poisson.rvs(lambdaBG)
        nSig=scipy.stats.poisson.rvs(lambdaSig)
        nObs=nBG+nSig
        relErrorBG=-1.
        while relErrorBG<0.:
            relErrorBG=scipy.stats.norm.rvs ( 0.3, 0.1 )
        estBG=scipy.stats.norm.rvs(nBG,math.sqrt ( nBG + ( relErrorBG*nBG )**2  ) )
      #  re=statistics.bayesianUpperLimit ( nObs, .0001, estBG, relErrorBG*estBG )
        print "nObs=",nObs
        print "nSig=",nSig
        print "nBG=",nBG
        print "estBG=",estBG,"+-",relErrorBG*estBG
        re=statistics.getUL ( nObs, estBG, relErrorBG*estBG )
        print "95% UL =",re


if __name__ == "__main__":
    unittest.main()
