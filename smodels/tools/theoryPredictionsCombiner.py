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
import scipy.stats as stats
import scipy.optimize as optimize

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
            lumi = tp.dataset.getLumi()
            nsig = (tp.xsection.value*lumi).asNumber()
            ret = ret * tp.dataset.likelihood(mu*nsig,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
        return ret

    def chi2( self ):
        if not hasattr ( self, "likelihood" ):
            logger.error ( "call computeStatistics before calling chi2" )
            return
        return - 2 * np.log ( self.likelihood / self.lmax )

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
            ## posterior expected
            self.elmax = lsm
        else:
            self.likelihood = llhd
            self.lsm = lsm
        self.mu_hat, self.sigma_mu, self.lmax = self.findMuHat ( marginalize = marginalize, deltas_rel = deltas_rel, extended_output = True )

    def findMuHat( self,marginalize=False, deltas_rel=.2, allowNegativeSignals = False, expected = False, extended_output = False, nll=False ):
        """ find muhat and lmax.
        :param allowNegativeSignals: if true, then also allow for negative values
        :param expected: if true, compute expected prior (=lsm), if "posteriori"
                         compute posteriori expected
        :param extended_output: if true, return also sigma_mu, the estimate of the error of mu_hat, and lmax, the likelihood at mu_hat
        :param nll: if true, return negative log max likelihood instead of lmax
        :returns: mu_hat, i.e. the maximum likelihood estimate of mu
        """
        import scipy.optimize
        def fun ( mu, *myargs ):
            return - np.log ( self.getLikelihood ( mu, marginalize = myargs[0],
                        deltas_rel = myargs[1], expected = expected, 
                        useRelSigStrengths = True ) )
        o = scipy.optimize.minimize ( fun, 1., args = ( marginalize, deltas_rel ) )
        lmax = np.exp ( - o.fun )
        logger.debug ( "result", o )
        mu_hat = o.x[0]
        if not allowNegativeSignals and mu_hat < 0.:
            mu_hat = 0.
        ## the inverted hessian is a good approximation for the variance at the
        ## minimum
        sigma_mu = np.sqrt ( o.hess_inv[0][0] )
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

    def getLikelihood(self,mu=1.,marginalize=False,deltas_rel=.2,expected=False,
                      nll = False, useRelSigStrengths = True ):
        """
        Compute the likelihood at a given mu
        :param mu: signal strength
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
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
            llhd = llhd * tp.dataset.likelihood(mu*nsig,marginalize=marginalize,
                                             deltas_rel=deltas_rel,expected=expected)
        if nll:
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

    def getUpperLimitOnMu(self, marginalize=False, expected=False, deltas_rel=0.2):
        """ get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: if True, compute expected likelihood, else observed
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :returns: upper limit on signal strength multiplier mu
        """
        if not hasattr ( self, "mu_hat" ):
            self.computeStatistics( marginalize = marginalize, 
                    deltas_rel = deltas_rel, expected = False )
        #nll0 = self.getLikelihood ( self.mu_hat, marginalize = marginalize,
        #                            nll = True )
        nll0 = self.getLikelihood ( self.mu_hat, marginalize = marginalize,
                        expected = expected, nll = True, useRelSigStrengths = True )
        toys = 30000
        #print ( f"COMB nll0 {nll0:.3f} mu_hat {self.mu_hat:.3f} sigma_mu {self.sigma_mu:.3f}" )
        ## a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        mu_hatA,_,nll0A = self.findMuHat ( marginalize = False, deltas_rel = .2,
                           expected = "posteriori", nll=True, extended_output= True )
        #print ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        #return 1.

        def root_func(mu):
            nll = self.getLikelihood( mu, marginalize=marginalize, nll=True,
                     expected= expected, useRelSigStrengths = True )
            nllA = self.getLikelihood( mu, marginalize=marginalize, 
                     expected="posteriori", nll=True, useRelSigStrengths = True )
            qmu =  2*( nll - nll0 )
            if qmu<0.: qmu=0.
            sqmu = np.sqrt (qmu)
            qA =  2*( nllA - nll0A )
            if qA<0.:
                qA = 0.
            sqA = np.sqrt(qA)
            CLsb = 1. - stats.multivariate_normal.cdf(sqmu)
            CLb = 0.
            if qA >= qmu:
                CLb =  stats.multivariate_normal.cdf(sqA - sqmu)
            else:
                if qA == 0.:
                    CLsb = 1.
                    CLb  = 1.
                else:
                    CLsb = 1. - stats.multivariate_normal.cdf( (qmu + qA)/(2*sqA) )
                    CLb = 1. - stats.multivariate_normal.cdf( (qmu - qA)/(2*sqA) )
            CLs = 0.
            if CLb > 0.:
                CLs = CLsb / CLb
            cl = .95
            root = CLs - 1. + cl
            return root

        a0 = root_func(0.)
        if a0 < 0. and marginalize:
            if toys < 20000:
                return self.getUpperLimitOnMu ( mu, marginalize, 4*toys,
                     expected = expected, trylasttime=False )
            else:
                if not trylasttime:
                    return self.getUpperLimitOnMu( mu, False, toys,
                                          expected = expected, trylasttime=True )
                # it ends here
                return None
        a,b=1.5*self.mu_hat,2.5*self.mu_hat+2*self.sigma_mu
        ctr=0
        while True:
            while ( np.sign ( root_func(a)* root_func(b) ) > -.5 ):
                b=1.7*b  ## widen bracket FIXME make a linear extrapolation!
                a=a-(b-a)*.4 ## widen bracket
                if a < 0.: a=0.
                ctr+=1
                if ctr>20: ## but stop after 20 trials
                    if toys > 20000 or not marginalize:
                        logger.error("cannot find brent bracket after 20 trials. f(a=%s)=%s,f(b=%s)=%s, mu_hat=%.2f, sigma_mu=%.2f" % ( a, root_func(a),b,root_func(b), self.mu_hat, self.sigma_mu ) )
                        logger.error("nll0=%s" % ( nll0 ) )
                        if marginalize == True:
                            return self.getUpperLimitOnMu( mu, marginalize = False,
                                                  toys = toys, expected = expected,
                                                  trylasttime = True )
                        return float("inf") ## better choice than None
                    else:
                        logger.debug("cannot find brent bracket after 20 trials. but very low number of toys")
                        return self.getUpperLimitOnMu ( mu, marginalize, 4*toys,
                                              expected=expected, trylasttime = False )
            try:
                # root_func bei 0 sollte positiv sein, bei inf = -.05
                mu_lim = optimize.brentq ( root_func, a, b, rtol=1e-03, xtol=1e-06 )
                # print ( f"COMBO mu_lim {mu_lim:.3f} expected {expected}" )
                return mu_lim 
            except ValueError as e: ## it could still be that the signs arent opposite
                # in that case, try again
                pass

        ## and so on and so on
        return None

    def getUpperLimitOnMuOld(self, marginalize=False, expected=False, deltas_rel=0.2):
        """ get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: if True, compute expected likelihood, else observed
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :returns: upper limit on signal strength multiplier mu
        """
        if not hasattr ( self, "mu_hat" ):
            self.computeStatistics( marginalize = marginalize, 
                    deltas_rel = deltas_rel, expected = False )
        ## lower and upper bounds on mu that we scan
        ## +/- 3 sigma should cover it up to < 1 permille
        lower, upper = self.mu_hat - 4*self.sigma_mu, self.mu_hat + 4*self.sigma_mu

        ## for one-sided cases we go up to 4.5 sigma, has similar p value.
        if lower < 0.:
            lower = 0.
            upper = 4.5 * self.sigma_mu
        if expected:
            lower,upper = 0., self.sigma_mu * 4.5
        llhds={}
        totllhd=0.
        nbins = 50.
        delta = upper / (nbins-1)
        for mu in np.arange ( lower, upper*1.0001, delta):
            llhd = self.getLikelihood ( mu, marginalize=marginalize, 
                                        deltas_rel = deltas_rel, expected=expected )
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

    def getUpperLimit(self, marginalize=False, expected=False, deltas_rel=0.2):
        """ get upper limit on *fiducial* cross sections
            which CLs = 0.95
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: if True, compute expected likelihood, else observed
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :returns: upper limit on *fiducial* cross sections (list)
        """
        clmu = self.getUpperLimitOnMu ( marginalize, expected, deltas_rel )
        ret = []
        for tp in self.theoryPredictions:
            ret.append ( tp.xsection.value * clmu )
        return ret

    def getRValue(self, marginalize=False, expected=False, deltas_rel=0.2):
        """ obtain r-value, i.e. predicted_xsec / ul
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: if True, compute expected likelihood, else observed
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        :returns: r-value
        """
        clmu = self.getUpperLimitOnMu ( marginalize, expected, deltas_rel )
        #ret = []
        xsecs = []
        for tp in self.theoryPredictions:
            xsecs.append ( float (tp.xsection.value * tp.dataset.getLumi() ) )
        #return ret
        #print ( "clmu", clmu )
        #print ( "xsecs", xsecs )
        return sum(xsecs) / clmu
