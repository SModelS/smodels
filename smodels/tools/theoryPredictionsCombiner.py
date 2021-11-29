#!/usr/bin/env python3

"""
.. module:: theoryPredictionsCombiner
   :synopsis: a module that deals with combining signal regions from
              different analyses, offering the same API as simplified likelihoods.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import math

class TheoryPredictionsCombiner():
    """
    Facility used to combine theory predictions from different analyes.
    """
    def __init__ (self, theoryPredictions : list ):
        self.theoryPredictions = theoryPredictions


    def lsm( self, marginalize = False, deltas_rel = .2, expected = False ):
        """ compute the likelihood at mu = 0.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :param expected: if True, compute expected likelihood, else observed
        """
        llhd = 1.
        for tp in self.theoryPredictions:
            llhd = llhd * tp.dataset.likelihood(0.,marginalize=marginalize,deltas_rel=deltas_rel,expected=expected)
        return llhd

    def likelihood ( self, mu = 1., marginalize = False, deltas_rel=.2, 
                     expected = False ):
        """ compute the likelihood at signal strength mu
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :param expected: if True, compute expected likelihood, else observed
        """
        ret = 1.
        for tp in self.theoryPredictions:
            lumi = self.dataset.getLumi()
            nsig = (self.xsection.value*lumi).asNumber()
            ret = ret * tp.dataset.likelihood(mu*nsig,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
        return ret

    def computeStatistics(self,marginalize=False,deltas_rel=0.2, expected=False ):
        """
        Compute the likelihoods, chi2, and upper limit for this combination
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :param expected: if True, compute expected likelihood, else observed
        """
        if expected == "both":
            for e in [ False, True ]:
                self.computeStatistics ( marginalize, deltas_rel, e )
            return
        llhd, lsm = 1., 1.
        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            llhd = llhd * tp.dataset.likelihood(nsig,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
            lsm = lsm * tp.dataset.likelihood(0.,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
        if expected:
            self.elikelihood = llhd
            self.elsm = lsm
            self.elmax = lsm
        else:
            self.likelihood = llhd
            self.lsm = lsm
        self.computeLMax ( marginalize = marginalize, deltas_rel = deltas_rel )

    def computeLMax(self,marginalize=False, deltas_rel=.2 ):
        """ find lmax. """
        import scipy.optimize
        def fun ( mu, *myargs ):
            return - math.log ( self.getLikelihood ( mu, marginalize = myargs[0],
                                                     deltas_rel = myargs[1] ) )
        o = scipy.optimize.minimize ( fun, 1., args = ( marginalize, deltas_rel ) )
        self.lmax = math.exp ( - o.fun )
        self.muhat = o.x[0]

    def getLikelihood(self,mu=1.,marginalize=False,deltas_rel=.2,expected=False):
        """
        Compute the likelihood at a given mu
        :param mu: signal strength
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :param expected: if True, compute expected likelihood, else observed
        """
        llhd = 1.
        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            llhd = llhd * tp.dataset.likelihood(mu*nsig,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
        return llhd

    def getUpperLimit(self, expected=False, deltas_rel=0.2):
        """ not yet implemented. get upper limit, i.e. value for mu for
            which CLs = 0.95
        """
        #if not hasattr ( self, "lsm" ):
        #    self.computeStatistics ( False, deltas_rel, False )
        lower,upper=self.muhat/5.,self.muhat*5.
        if expected:
            lower,upper = 0., self.muhat*3.
        llhds={}
        totllhd=0.
        for mu in np.arange ( lower, upper, lower ):
            llhd = self.getLikelihood ( mu, 

        return None

    def getRValue ( self, expected = False ):
        """ not yet implemented, get r-value """
        upperLimit = self.getUpperLimit(expected)
        if upperLimit is None or upperLimit.asNumber(fb)==0.:
            return None
        return None

    def getPValue ( self ):
        """ is this needed? """
        return None
