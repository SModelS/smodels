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
import scipy.optimize as optimize
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError


class TheoryPredictionsCombiner(object):
    """
    Facility used to combine theory predictions from different analyes.
    """

    def __init__(self, theoryPredictions: list,
                 slhafile=None, marginalize=False, deltas_rel=None):
        """ constructor.
        :param theoryPredictions: the List of theory predictions
        :param slhafile: optionally, the slhafile can be given, for debugging
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        if len(theoryPredictions) == 0:
            raise SModelSError("asking to combine zero theory predictions")
        self.theoryPredictions = theoryPredictions
        self.slhafile = slhafile
        self.marginalize = marginalize
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default
            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedLlhds = {False: {}, True: {}, "posteriori": {}}

    def __str__(self):
        ret = "%s:%s" % (self.analysisId(), self.totalXsection())
        return ret

    def singleDecorator(function):
        """
        If the combiner holds a single theoryPrediction, calls
        the theoryPrediction function.
        """

        def wrapper(self, *args, **kwargs):
            # Check if obj contains a single prediction
            if len(self.theoryPredictions) == 1:
                tp = self.theoryPredictions[0]
                # Check if the theoryPrediction object has the function
                if hasattr(tp, function.__name__):
                    newF = getattr(tp, function.__name__)
                    # Check if the attribute is callable:
                    if callable(newF):
                        # Return value from theoryPrediction.function call
                        return newF(*args, **kwargs)

            # If anything failed, call method from combiner:
            return function(self, *args, **kwargs)

        return wrapper

    @singleDecorator
    def dataId(self):
        """
        Return a string with the IDs of all the datasets used in the combination.
        """
        ret = ','.join([tp.dataset.getID() for tp in self.theoryPredictions])

        return ret

    @singleDecorator
    def analysisId(self):
        """
        Return a string with the IDs of all the experimental results used in the combination.
        """

        ret = ','.join(sorted([tp.expResult.globalInfo.id for tp in self.theoryPredictions]))

        return ret

    @singleDecorator
    def dataType(self, short=False):
        """
        Return the type of dataset (=combined)
        :param: short, if True, return abbreviation (comb)
        """
        if short:
            return 'comb'
        else:
            return 'combined'

    @singleDecorator
    def getUpperLimit(self, expected=False,
                      trylasttime=False):
        """ get upper limit on *fiducial* cross sections
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :returns: upper limit on *fiducial* cross sections (list)
        """
        clmu = self.getUpperLimitOnMu(expected=expected)
        ret = []
        for tp in self.theoryPredictions:
            ret.append(tp.xsection.value * clmu)
        return ret

    @singleDecorator
    def getRValue(self, expected=False):
        """ obtain r-value, i.e. predicted_xsec / ul
        :param expected: if True, compute expected likelihood, else observed
        :returns: r-value
        """
        if "r" in self.cachedObjs[expected]:
            return self.cachedObjs[expected]["r"]
        clmu = self.getUpperLimitOnMu(expected=expected)
        if clmu == 0.:
            # has no meaning, has it?
            return None
        r = 1.0 / clmu
        self.cachedObjs[expected]["r"] = r
        return r

    @singleDecorator
    def lsm(self, expected=False):
        """ compute the likelihood at mu = 0.
        :param expected: if True, compute expected likelihood, else observed.
                         If "posteriori", compute posterior expected.
        """
        llhd = 1.
        for tp in self.theoryPredictions:
            llhd = llhd * tp.likelihood(0., expected=expected)
        return llhd

    @singleDecorator
    def lmax(self, expected=False):
        if "lmax" not in self.cachedObjs[expected]:
            self.computeStatistics(expected=expected)
        if "lmax" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"] = None
        return self.cachedObjs[expected]["lmax"]

    @singleDecorator
    def muhat(self, expected=False):
        if "muhat" not in self.cachedObjs[expected]:
            self.computeStatistics(expected=expected)
        if "muhat" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["muhat"] = None
        return self.cachedObjs[expected]["muhat"]

    @singleDecorator
    def chi2(self, expected=False):
        if "llhd" not in self.cachedObjs[expected] or "lmax" not in self.cachedObjs[expected]:
            logger.error("call computeStatistics before calling chi2")
            return
        llhd = self.cachedObjs[expected]["llhd"]
        lmax = self.cachedObjs[expected]["lmax"]

        if llhd == 0.:
            return 2000.  # we cut off at > 1e-300 or so ;)
        return - 2 * np.log(llhd / lmax)

    @singleDecorator
    def likelihood(self, mu=1., expected=False,
                   nll=False, useCached=True):
        """
        Compute the likelihood at a given mu
        :param mu: signal strength
        :param expected: if True, compute expected likelihood, else observed
        :param nll: if True, return negative log likelihood, else likelihood
        :param useCached: if True, will return the cached value from theoryPrediction (if available)
        """
        try:
            mu = mu[0]  # some of these methods use arrays with a single element
        except:
            pass
        if mu in self.cachedLlhds[expected]:
            llhd = self.cachedLlhds[expected][mu]
            if nll:
                if llhd == 0.:  # cut off nll at 700 (1e-300)
                    return 700.
                return - np.log(llhd)
            return llhd

        llhd = 1.

        for tp in self.theoryPredictions:
            # Set the theory marginalize attribute to the combiner value:
            tp_marginalize = tp.marginalize
            tp.marginalize = self.marginalize
            llhd = llhd * tp.likelihood(mu, expected=expected, useCached=useCached)
            # Restore marginalize setting:
            tp.marginalize = tp_marginalize
        self.cachedLlhds[expected][mu] = llhd
        if nll:
            if llhd == 0.:  # cut off nll at 700 (1e-300)
                return 700.
            return - np.log(llhd)
        return llhd

    @singleDecorator
    def computeStatistics(self, expected=False):
        """
        Compute the likelihoods, chi2, and upper limit for this combination
        :param expected: if True, compute expected likelihood, else observed
        """
        if self.theoryPredictions is None:
            return
        if expected == "both":
            for e in [False, True]:
                self.computeStatistics(e)
            return
        llhd, lsm = 1., 1.
        for tp in self.theoryPredictions:
            tp_marginalize = tp.marginalize
            tp.marginalize = self.marginalize
            llhd = llhd * tp.likelihood(1.0, expected=expected,
                                        useCached=True)

            lsm = lsm * tp.likelihood(0.0, expected=expected,
                                      useCached=True)
            # Restore marginalize setting:
            tp.marginalize = tp_marginalize

        self.cachedObjs[expected]["llhd"] = llhd
        self.cachedObjs[expected]["lsm"] = lsm
        if expected:
            self.cachedObjs[expected]["lmax"] = lsm
        self.mu_hat, self.sigma_mu, lmax = self.findMuHat(extended_output=True)
        self.cachedObjs[expected]["lmax"] = lmax
        self.cachedObjs[expected]["muhat"] = self.mu_hat

    @singleDecorator
    def totalXsection(self):
        ret = 0.*fb
        if self.theoryPredictions is not None:
            for tp in self.theoryPredictions:
                ret += tp.xsection.value
        return ret

    @singleDecorator
    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """
        conditions = [tp.getmaxCondition() for tp in self.theoryPredictions]
        return max(conditions)

    def findMuHat(self, allowNegativeSignals=False, expected=False,
                  extended_output=False, nll=False):
        """ find muhat and lmax.
        :param allowNegativeSignals: if true, then also allow for negative values
        :param expected: if true, compute expected prior (=lsm), if "posteriori"
                         compute posteriori expected
        :param extended_output: if true, return also sigma_mu, the estimate of the error of mu_hat, and lmax, the likelihood at mu_hat
        :param nll: if true, return negative log max likelihood instead of lmax
        :returns: mu_hat, i.e. the maximum likelihood estimate of mu
        """
        import scipy.optimize

        def fun(mu):
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            return self.likelihood(mu, expected=expected,
                                   nll=True, useCached=False)

        toTry = [1., 0., 3., 10.]
        if allowNegativeSignals:
            toTry = [1., 0., 3., -1., 10., -3.]
        for mu0 in toTry:
            # Minimize with a delta_mu = 1e-3*mu0 when computing the first derivative
            # (if delta_mu is too small, the derivative might give zero and terminate the minimization)
            o = scipy.optimize.minimize(fun, mu0)
            logger.debug(f"result for {mu0}: %s" % o)
            if not o.success:
                logger.debug(
                    f"combiner.findMuHat did not terminate successfully: {o.message} mu_hat={o.x} hess={o.hess_inv} slhafile={self.slhafile}")
            # the inverted hessian is a good approximation for the variance at the
            # minimum
            hessian = o.hess_inv[0][0]
            if hessian > 0.:  # found a maximum. all good.
                break
            # the hessian is negative meaning we found a maximum, not a minimum
            logger.debug(
                f"combiner.findMuHat the hessian {hessian} is negative at mu_hat={o.x} in {self.slhafile} try again with different initialisation.")
        mu_hat = o.x[0]
        lmax = np.exp(- o.fun)
        if hessian < 0.:
            logger.error(
                "tried with several starting points to find maximum, always ended up in minimum. bailing out.")
            if extended_output:
                return None, None, None
            return None
        if not allowNegativeSignals and mu_hat < 0.:
            mu_hat = 0.  # fixme for this case we should reevaluate the hessian!
        sigma_mu = np.sqrt(hessian)
        if nll:
            lmax = - np.log(lmax)
        if extended_output:
            return mu_hat, sigma_mu, lmax

        return mu_hat

    def getUpperLimitOnMu(self, expected=False):
        """ get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :returns: upper limit on signal strength multiplier mu
        """
        if "UL" in self.cachedObjs[expected]:
            return self.cachedObjs[expected]["UL"]
        # if not hasattr ( self, "mu_hat" ):
        #    self.computeStatistics( expected = False )
        mu_hat, sigma_mu, lmax = self.findMuHat(allowNegativeSignals=True,
                                                extended_output=True)
        nll0 = self.likelihood(mu_hat, expected=expected, nll=True)
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        mu_hatA, _, nll0A = self.findMuHat(expected="posteriori", nll=True,
                                           extended_output=True)
        # print ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.

        def clsRoot(mu):
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll = self.likelihood(mu, nll=True, expected=expected,
                                  useCached=False)
            nllA = self.likelihood(mu, expected="posteriori", nll=True,
                                   useCached=False)
            ret = rootFromNLLs(nllA, nll0A, nll, nll0)
            return ret

        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot)
        mu_lim = optimize.brentq(clsRoot, a, b, rtol=1e-03, xtol=1e-06)
        self.cachedObjs[expected]["UL"] = mu_lim
        return mu_lim
