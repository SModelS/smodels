#!/usr/bin/env python3

"""
.. module:: theoryPredictionsCombiner
   :synopsis: a module that deals with combining signal regions from
              different analyses, offering the same API as simplified likelihoods.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import numpy as np
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.tools.statistics import rootFromNLLs, determineBrentBracket
import scipy.stats as stats
import scipy.optimize as optimize

class TheoryPredictionsCombiner():
    """
    Facility used to combine theory predictions from different analyes.
    """
    def __init__ ( self, theoryPredictions : list, slhafile = None,
                   toys = 30000, marginalize = False, deltas_rel = None ):
        """ constructor.
        :param theoryPredictions: the List of theory predictions
        :param slhafile: optionally, the slhafile can be given, for debugging
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        self.theoryPredictions = theoryPredictions
        self.slhafile = slhafile
        self.toys = toys
        self.marginalize = marginalize
        if deltas_rel == None:
            from smodels.tools.runtime import _deltas_rel_default
            deltas_rel = _deltas_rel_default
        self.deltas_rel  = deltas_rel
        self.cachedObjs = { False: {}, True: {}, "posteriori": {} }

    def lsm( self, expected = False ):
        """ compute the likelihood at mu = 0.
        :param expected: if True, compute expected likelihood, else observed.
                         If "posteriori", compute posterior expected.
        """
        llhd = 1.
        for tp in self.theoryPredictions:
            llhd = llhd * tp.dataset.likelihood(0.,marginalize=self.marginalize,
                    deltas_rel=self.deltas_rel,expected=expected )
        return llhd

    def likelihood ( self, mu = 1., expected = False ):
        """ compute the likelihood at signal strength mu
        :param expected: if True, compute expected likelihood, else observed
        """
        ret = 1.
        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            ret = ret * tp.dataset.likelihood(mu*nsig,marginalize=self.marginalize,
                                        deltas_rel=self.deltas_rel,expected=expected)
        return ret

    def likelihood( self, expected=False ):
        return self.cachedObjs[expected]["llhd"]

    def lmax( self, expected=False ):
        return self.cachedObjs[expected]["lmax"]

    def chi2( self, expected=False ):
        if not "llhd" in self.cachedObjs[expected] or not "lmax" in self.cachedObjs[expected]:
            logger.error ( "call computeStatistics before calling chi2" )
            return
        llhd = self.cachedObjs[expected]["llhd"]
        lmax = self.cachedObjs[expected]["lmax"]
        if llhd == 0.:
            return 2000. ## we cut off at > 1e-300 or so ;)
        return - 2 * np.log ( llhd / lmax )

    def totalXsection ( self ):
        ret = 0.*fb
        if self.theoryPredictions != None:
            for tp in self.theoryPredictions:
                ret += tp.xsection.value
        return ret

    def getmaxCondition ( self ):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """
        conditions = [ tp.getmaxCondition() for tp in self.theoryPredictions ]
        return max(conditions)


    def computeStatistics(self, expected=False ):
        """
        Compute the likelihoods, chi2, and upper limit for this combination
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :param expected: if True, compute expected likelihood, else observed
        """
        if self.theoryPredictions == None:
            return
        if expected == "both":
            for e in [ False, True ]:
                self.computeStatistics ( e )
            return
        llhd, lsm = 1., 1.
        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            llhd = llhd * tp.dataset.likelihood(nsig,marginalize=self.marginalize,
                                       deltas_rel=self.deltas_rel,expected=expected)
            lsm = lsm * tp.dataset.likelihood(0.,marginalize=self.marginalize,
                                       deltas_rel=self.deltas_rel,expected=expected)
        self.cachedObjs[expected]["llhd"]=llhd
        self.cachedObjs[expected]["lsm"]=lsm
        if expected:
            self.cachedObjs[expected]["lmax"]=lsm
        self.mu_hat, self.sigma_mu, lmax = self.findMuHat ( allowNegativeSignals = False, extended_output = True )
        self.cachedObjs[expected]["lmax"]=lmax

    def findMuHat( self, allowNegativeSignals = False, expected = False, 
                   extended_output = False, nll=False ):
        """ find muhat and lmax.
        :param allowNegativeSignals: if true, then also allow for negative values
        :param expected: if true, compute expected prior (=lsm), if "posteriori"
                         compute posteriori expected
        :param extended_output: if true, return also sigma_mu, the estimate of the error of mu_hat, and lmax, the likelihood at mu_hat
        :param nll: if true, return negative log max likelihood instead of lmax
        :returns: mu_hat, i.e. the maximum likelihood estimate of mu
        """
        import scipy.optimize
        def fun ( mu ):
            return - np.log ( self.getLikelihood ( mu, expected = expected, 
                        useRelSigStrengths = True ) )
        toTry = [ 1., 0., 3., 10. ]
        if allowNegativeSignals:
            toTry = [ 1., 0., 3., -1., 10., -3. ]
        for mu0 in toTry:
            o = scipy.optimize.minimize ( fun, mu0 )
            logger.debug ( f"result for {mu0}:", o )
            if not o.success:
                logger.debug ( f"combiner.findMuHat did not terminate successfully: {o.message} mu_hat={o.x} hess={o.hess_inv} slhafile={self.slhafile}" )
            ## the inverted hessian is a good approximation for the variance at the
            ## minimum
            hessian = o.hess_inv[0][0]
            if hessian>0.: # found a maximum. all good.
                break
            # the hessian is negative meaning we found a maximum, not a minimum
            logger.debug ( f"combiner.findMuHat the hessian {hessian} is negative at mu_hat={o.x} in {self.slhafile} try again with different initialisation." )
        mu_hat = o.x[0]
        lmax = np.exp ( - o.fun )
        if hessian < 0.:
            logger.error ( "tried with several starting points to find maximum, always ended up in minimum. bailing out." )
            if extended_output:
                return None, None, None
            return None
        if not allowNegativeSignals and mu_hat < 0.:
            mu_hat = 0. # fixme for this case we should reevaluate the hessian!
        sigma_mu = np.sqrt ( hessian )
        if nll:
            lmax = - np.log(lmax)
        if extended_output:
            return mu_hat, sigma_mu, lmax
        return mu_hat
    
    def getNSigtot ( self ):
        """ get the total number of signal events, summed over
            all theory predictions. needed for normalizations. """
        S = 0.
        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            S += nsig
        return S

    def getLikelihood(self,mu=1.,expected=False,
                      nll = False, useRelSigStrengths = True ):
        """
        Compute the likelihood at a given mu
        :param mu: signal strength
        :param expected: if True, compute expected likelihood, else observed
        :param nll: if True, return negative log likelihood, else likelihood
        :param useRelSigStrengths: if True, get the likelihood not on signal strength mu, but on relative signal strength, total_nsig = 1 for mu = 1
        """
        llhd = 1.
        S = 1.
        if useRelSigStrengths:
            S = self.getNSigtot()

        for tp in self.theoryPredictions:
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            if useRelSigStrengths:
                nsig = nsig / S
            llhd = llhd * tp.dataset.likelihood(mu*nsig,marginalize=self.marginalize,
                                        deltas_rel=self.deltas_rel,expected=expected)
        if nll:
            if llhd == 0.: ## cut off nll at 700 (1e-300)
                return 700.
            return - np.log ( llhd )
        return llhd

    def signals ( self, mu ):
        """
        Returns the number of expected signal events, for all datasets,
        given total signal strength mu.
        :param mu: Total number of signal events summed over all datasets.
        """
        ret = []
        for tp in self.theoryPredictions:
            n = float ( mu * tp.dataset.getLumi() * tp.xsection.value )
            ret.append ( n )
        return ret

    def getUpperLimitOnMu(self, toys = None, expected = False ):
        """ get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :param toys: specify number of toys. Use default is none
        :returns: upper limit on signal strength multiplier mu
        """
        if toys==None:
            toys=self.toys
        #if not hasattr ( self, "mu_hat" ):
        #    self.computeStatistics( expected = False )
        mu_hat, sigma_mu, lmax = self.findMuHat ( allowNegativeSignals = True, extended_output = True )
        nll0 = self.getLikelihood ( mu_hat, expected = expected, nll = True, 
                useRelSigStrengths = True )
        ## a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        mu_hatA,_,nll0A = self.findMuHat ( expected = "posteriori", nll=True, 
                extended_output= True )
        #print ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        #return 1.

        def clsRoot (mu):
            # at - infinity this should be .95,
            # at + infinity it should -.05
            nll = self.getLikelihood( mu, nll=True,
                     expected= expected, useRelSigStrengths = True )
            nllA = self.getLikelihood( mu,  expected="posteriori", nll=True, 
                    useRelSigStrengths = True )
            return rootFromNLLs ( nllA, nll0A, nll, nll0 )

        a,b = determineBrentBracket( mu_hat, sigma_mu, clsRoot )
        mu_lim = optimize.brentq ( clsRoot, a, b, rtol=1e-03, xtol=1e-06 )
        return mu_lim 

    def getUpperLimitOnMuOld(self, expected=False ):
        """ get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :returns: upper limit on signal strength multiplier mu
        """
        if not hasattr ( self, "mu_hat" ):
            self.computeStatistics( marginalize = marginalize, 
                    deltas_rel = self.deltas_rel, expected = False )
        ## lower and upper bounds on mu that we scan
        ## +/- 3 sigma should cover it up to < 1 permille
        lower, upper = self.mu_hat - 4*sigma_mu, mu_hat + 4*sigma_mu

        ## for one-sided cases we go up to 4.5 sigma, has similar p value.
        if lower < 0.:
            lower = 0.
            upper = 4.5 * sigma_mu
        if expected:
            lower,upper = 0., sigma_mu * 4.5
        llhds={}
        totllhd=0.
        nbins = 50.
        delta = upper / (nbins-1)
        for mu in np.arange ( lower, upper*1.0001, delta):
            llhd = self.getLikelihood ( mu, expected=expected )
            llhds[mu]=llhd
            totllhd+=llhd
        ## last likelihood should contribute less than 1 permille
        if llhd > totllhd * 1e-3:
            logger.warning ( f"when discretizing likelihood, last bin contributes too much {llhd/totllhd}" )

        cum, lastcum, lastmu = 0., 0., 0.
        alpha = .05
        significance = 1. - alpha
        for mu_,cdf_ in llhds.items():
            cum += cdf_/totllhd
            llhds[mu_]=cdf_/totllhd
            if cum > significance:
                k = ( cum - lastcum ) / ( mu_ - lastmu )
                d = cum - k * mu_
                mu_ret = ( significance - d ) / k
                return mu_ret
            lastcum = cum
            lastmu = mu_

    def getUpperLimit( self, toys =None, expected=False, 
                       trylasttime = False ):
        """ get upper limit on *fiducial* cross sections
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :returns: upper limit on *fiducial* cross sections (list)
        """
        clmu = self.getUpperLimitOnMu ( expected = expected )
        ret = []
        for tp in self.theoryPredictions:
            ret.append ( tp.xsection.value * clmu )
        return ret

    def getRValue(self, expected=False ):
        """ obtain r-value, i.e. predicted_xsec / ul
        :param expected: if True, compute expected likelihood, else observed
        :returns: r-value
        """
        clmu = self.getUpperLimitOnMu ( expected = expected )
        if clmu == 0.:
            ## has no meaning, has it?
            return None
        #ret = []
        xsecs = []
        for tp in self.theoryPredictions:
            xsecs.append ( float (tp.xsection.value * tp.dataset.getLumi() ) )
        #return ret
        #print ( "clmu", clmu )
        #print ( "xsecs", xsecs )
        return sum(xsecs) / clmu
