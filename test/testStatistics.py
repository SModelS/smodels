#!/usr/bin/env python

"""
.. module:: testStatistics
   :synopsis: Tests the statistics module.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools import statistics
from smodels.tools.physicsUnits import fb
import scipy.stats
import math
import random

class StatisticsTest(unittest.TestCase):
    def testUpperLimit(self):
        re = statistics.upperLimit ( 100., 100., 0., 20./fb   )
        self.assertAlmostEqual ( re.asNumber ( fb ), 1.06, 1 )
    """
    def testCLInterval(self):
        re = statistics.computeCLInterval ( 100., 100., 1. )
        self.assertAlmostEqual ( re, 18.0792727821 )

    def testBayesianLimit(self):
        re = statistics.bayesianUpperLimit ( 100, 0., 100., 0. )
        self.assertAlmostEqual ( re, 21.4256127382 )

    def testUL(self):
        re = statistics.getUL ( 100, 100., 0. )

    def westCoverage(self):
        coverage=[]
        for i in range(100):
            lambdaBG=random.uniform(10,50)
            lambdaSig=random.uniform(5,30)
            nBG=scipy.stats.poisson.rvs(lambdaBG)
            nSig=scipy.stats.poisson.rvs(lambdaSig)
            nObs=nBG+nSig
            relErrorBG=-1.
            estBG=-1.
            while relErrorBG<0.:
                relErrorBG=scipy.stats.norm.rvs ( 0.2, 0.15 )
            relErrorBG=0.0001
            while estBG<0.:
                estBG=scipy.stats.norm.rvs( lambdaBG,math.sqrt ( lambdaBG + ( relErrorBG*lambdaBG )**2  ) )
                # estBG=scipy.stats.norm.rvs( lambdaBG,math.sqrt ( lambdaBG + ( relErrorBG*lambdaBG )**2  ) )
            #print "lambdaBG=",lambdaBG
            #print "nObs=",nObs
            #print "nSig=",nSig
            #print "nBG=",nBG
            #print "estBG=",estBG,"+-",relErrorBG*estBG
            try:
                re=statistics.getUL ( nObs, estBG, relErrorBG*estBG )
            except Exception,e:
                continue
            # re=statistics.bayesianUpperLimit ( nObs, .0001, estBG, relErrorBG*estBG )
            if re==0.0:
                continue
            #print "95% UL =",re
            #print "---------------"
            coverage.append ( re>nSig )
        print "coverage=",sum(coverage),"/",len(coverage)
        """



if __name__ == "__main__":
    unittest.main()
