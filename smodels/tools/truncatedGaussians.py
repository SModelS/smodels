#!/usr/bin/env python3

"""
.. module:: truncatedGaussian
   :synopsis: a module that contains the code that computes an approximate
   Gaussian likelihood from an expected an observer upper limit. See
   https://arxiv.org/abs/1202.3415.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
"""

__all__ = [ "TruncatedGaussians" ]

from scipy import stats, optimize
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from scipy.special import erf
import numpy as np
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from typing import Text, Optional, Union
from smodels.tools.basicStats import deltaChi2FromLlhd

class TruncatedGaussians:
    """ likelihood computer based on the trunacated Gaussian approximation, see
         arXiv:1202.3415 """

    def __init__  ( self, upperLimit, expectedUpperLimit, predicted_yield, 
                    corr : Optional[float] = 0.6, cl=.95, lumi = None ):
        """
        :param upperLimit: observed upper limit, as a yield or on xsec
        :param expectedUpperLimit: expected upper limit, also as a yield or on xsec
        :param predicted_yield: the predicted signal yield, unitless or [fb]
        :param corr: correction factor:
           ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
           When comparing with likelihoods constructed from efficiency maps,
           a factor of corr = 0.6 has been found to result in the best approximations.
        :param cl: confidence level
        :param lumi: if not None, and if the limits are in [fb], then use
                     it to translate limits on xsecs into limits on yields 
                     internally
        """
        if type(lumi) != type(None) and type(upperLimit) == type(fb):
            upperLimit = float ( upperLimit * lumi )
            expectedUpperLimit = float ( expectedUpperLimit * lumi )
            predicted_yield = float ( predicted_yield * lumi ) # the xsec
        if corr > 0.0 and upperLimit > expectedUpperLimit:
            expectedUpperLimit = expectedUpperLimit / (
                1.0 - corr * ((upperLimit - expectedUpperLimit) / (upperLimit + expectedUpperLimit))
            )
        self.upperLimit = upperLimit
        self.expectedUpperLimit = expectedUpperLimit
        self.predicted_yield = predicted_yield
        self.corr = corr
        self.sigma_y = self.getSigmaY()  # the expected scale, eq 3.24 in arXiv:1202.3415
        self.denominator = np.sqrt(2.0) * self.sigma_y
        self.cl = cl

    def likelihood ( self, mu : Union[float,None], nll : Optional[bool]=False, 
            allowNegativeMuhat : Optional[bool] = True,
            corr : Optional[float] = 0.6 ) -> float:
        """ return the likelihood, as a function of mu
        :param mu: number of signal events, if None then mu = muhat
        :param nll: if True, return negative log likelihood
        :param allowNegativeMuhat: if True, then allow muhat to become negative,
               else demand that muhat >= 0. In the presence of underfluctuations
               in the data, setting this to True results in more realistic
               approximate likelihoods.

        :returns: likelihood (float), muhat, and sigma_mu
        """
        sllhd = "llhd"
        if nll:
            sllhd = "nll"
        nsig = mu
        muhat, sigma_mu = float("inf"), float("inf")
        if mu != None:
            nsig = mu * self.predicted_yield
        dsig = self.likelihoodOfNSig ( nsig, nll=nll,
                allowNegativeMuhat = allowNegativeMuhat, corr = corr )
        if self.predicted_yield > 0.:
             muhat, sigma_mu =  dsig["yhat"]/self.predicted_yield,\
                dsig["sigma_y"] / self.predicted_yield

        ret = { sllhd: dsig[sllhd], "muhat": muhat, "sigma_mu": sigma_mu }
        return ret

    def likelihoodOfNSig ( self, nsig : Union[float,None], nll : Optional[bool]=False, 
            allowNegativeMuhat : Optional[bool] = True,
            corr : Optional[float] = 0.6 ) -> float:
        """ return the likelihood, as a function of nsig
        :param nsig: number of signal events, if None then nsig = muhat * predicted_yioelds
        :param nll: if True, return negative log likelihood
        :param allowNegativeMuhat: if True, then allow muhat to become negative,
               else demand that muhat >= 0. In the presence of underfluctuations
               in the data, setting this to True results in more realistic
               approximate likelihoods.
        :param corr: correction factor:
                     ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
                     When comparing with likelihoods constructed from efficiency maps,
                     a factor of corr = 0.6 has been found to result in the best approximations.

        :returns: likelihood (float), yhat, and sigma_y
        """
        sllhd = "llhd"
        if nll:
            sllhd = "nll"
    
        if self.upperLimit < self.expectedUpperLimit:
            ## underfluctuation. muhat = 0.
            if allowNegativeMuhat:
                xa = -self.expectedUpperLimit
                xb = 1
                yhat = self.find_neg_yhat( xa, xb )
                self.llhd_ = self.llhd(nsig, yhat, nll = False )
                ret = self.llhd_
                if nll:
                    ret = -math.log(ret)
                return { sllhd: ret, "yhat": yhat, "sigma_y": self.sigma_y }
            else:
                self.llhd_ = self.llhd(nsig, 0.0, nll = False )
                ret = self.llhd_
                if nll:
                    ret = -math.log(ret)
                return { sllhd: ret, "yhat": 0.0, "sigma_y": self.sigma_y }

        yhat = self.findYhat()
        self.llhd_ = self.llhd(nsig, yhat, nll = False )
        ret = self.llhd_
        if nll:
            ret = -math.log(ret)
        return { sllhd: ret, "yhat": yhat, "sigma_y": self.sigma_y }

    def findYhat ( self ):
        """ find the signal yields that maximize the likelihood """
        fA = self.root_func(0.0)
        fB = self.root_func(max(self.upperLimit, self.expectedUpperLimit))
        if np.sign(fA * fB) > 0.0:
            ## the have the same sign
            logger.error("when computing likelihood: fA and fB have same sign")
            return None, None, None
        yhat = optimize.brentq(
            self.root_func, 0.0, max(self.upperLimit, self.expectedUpperLimit), 
            rtol=1e-03, xtol=1e-06)
        return yhat

    def getSigmaY(self, yhat=0.0 ):
        """get the standard deviation sigma, given
        an upper limit and a central value. assumes a truncated Gaussian likelihood"""
        # the expected scale, eq 3.24 in arXiv:1202.3415
        return ( self.expectedUpperLimit - yhat) / 1.96

    def root_func(self,x):  # we want the root of this one
        return (erf((self.upperLimit - x) / self.denominator) + erf(x / self.denominator)) / (
            1.0 + erf(x / self.denominator)) - self.cl

    def find_neg_yhat(self, xa, xb):
        c = 0
        while self.root_func(xa) * self.root_func(xb) > 0:
            xa = 2 * xa
            c += 1
            if c > 10:
                logger.error(
                    f"cannot find bracket for brent bracketing ul={upperLimit}, eul={expectedUpperLimit},xa={xa}, xb={xb}"
                )

        muhat = optimize.brentq(self.root_func, xa, xb, rtol=1e-03, xtol=1e-06)
        return muhat

    def llhd( self, nsig, muhat, nll):
        if nsig is None:
            nsig = muhat
        # need to account for the truncation!
        # first compute how many sigmas left of center is 0.
        Zprime = muhat / self.sigma_y
        # now compute the area of the truncated gaussian
        A = stats.norm.cdf(Zprime)
        if nll:
            return np.log(A) - stats.norm.logpdf(nsig, muhat, self.sigma_y)
        return float(stats.norm.pdf(nsig, muhat, self.sigma_y) / A)

    def chi2( self, likelihood = None ):
        """compute the chi2 value from a likelihood (convenience function).
        :param likelihood: supply likelihood, if None, use just calculcated llhd
        """
        if likelihood == None:
            if not hasattr ( self, "llhd_" ):
                raise SModelSError ( "asking for chi2 but no likelihood given" )
            likelihood = self.llhd_
        l0 = 2.0 * stats.norm.logpdf(0.0, 0.0, self.sigma_y)
        l = deltaChi2FromLlhd(likelihood)
        if l is None:
            return None

        return l + l0

    def rvsFromLimits( self, n=1 ):
        """
        Generates a sample of random variates, given expected and observed likelihoods.
        The likelihood is modelled as a truncated Gaussian.

        :param n: sample size

        :returns: sample of random variates
        """
        muhat = self.findMuhat()
        ret = []
        while len(ret) < n:
            tmp = stats.norm.rvs(muhat, self.sigma_y)
            if tmp > 0.0:
                ret.append(tmp)
        return ret
