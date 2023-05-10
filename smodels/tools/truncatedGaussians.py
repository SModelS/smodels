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
from typing import Text, Optional, Union, Dict
from smodels.tools.basicStats import deltaChi2FromLlhd

class TruncatedGaussians:
    """ likelihood computer based on the trunacated Gaussian approximation, see
         arXiv:1202.3415 """
    __slots__ = [ "upperLimitOnMu", "expectedUpperLimitOnMu", "corr",
                  "sigma_mu", "denominator", "cl" ]

    # the correction we introduced a long time ago was have exclusions
    # in mind only. for discovery mode we need to correct differently
    newCorrectionType = False

    def __init__  ( self, upperLimitOnMu : float, expectedUpperLimitOnMu : float,
                    corr : Optional[float] = 0.6, cl=.95 ):
        """
        :param upperLimitOnMu: observed upper limit on signal strength mu
        :param expectedUpperLimitOnMu: expected upper limit on signal strength mu
        :param corr: correction factor:
           ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
           When comparing with likelihoods constructed from efficiency maps,
           a factor of corr = 0.6 has been found to result in the best approximations.
        :param cl: confidence level
        """
        assert type(upperLimitOnMu) in [ float, np.float64, np.float32 ], f"the upper limits must be given as floats not {type(upperLimitOnMu)}, are you providing upper limits on xsecs?"
        # the old type of correction
        if not self.newCorrectionType and corr > 0.0 and upperLimitOnMu > expectedUpperLimitOnMu:
            f = 1.0 - corr * ((upperLimitOnMu - expectedUpperLimitOnMu) / (upperLimitOnMu + expectedUpperLimitOnMu))
            expectedUpperLimitOnMu = expectedUpperLimitOnMu / f
        self.upperLimitOnMu = upperLimitOnMu
        self.expectedUpperLimitOnMu = expectedUpperLimitOnMu
        self.corr = corr
        self.sigma_mu = self._getSigmaMu()  # the expected scale, eq 3.24 in arXiv:1202.3415
        self.denominator = np.sqrt(2.0) * self.sigma_mu
        self.cl = cl

    def likelihood ( self, mu : Union[float,None], return_nll : Optional[bool]=False,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            expected : Union[Text,bool] = False ) -> Union[None,float]:
        """ return the likelihood, as a function of mu
        :param mu: number of signal events, if None then mu = muhat
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: if True, then allow muhat to become negative,\
               else demand that muhat >= 0. In the presence of underfluctuations\
               in the data, setting this to True results in more realistic\
               approximate likelihoods.

        :returns: likelihood (float)
        """
        sllhd = "llhd"
        if return_nll:
            sllhd = "nll"
        muhat, sigma_mu = float("inf"), float("inf")
        dsig = self._likelihoodOfMu ( mu, return_nll=return_nll,
                allowNegativeSignals = allowNegativeSignals, corr = corr,
                expected = expected )
        ret = dsig[sllhd]
        return ret

    def lmax ( self, return_nll : Optional[bool]=False,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            expected : Union[bool,Text] = False ) -> Dict:
        """ return the likelihood, as a function of mu
        :param mu: number of signal events, if None then mu = muhat
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: if True, then allow muhat to become negative,
               else demand that muhat >= 0. In the presence of underfluctuations
               in the data, setting this to True results in more realistic
               approximate likelihoods.

        :returns: dictionary with likelihood (float), muhat, and sigma_mu
        """
        default = { "muhat": None, "sigma_mu": None, "lmax": None }
        sllhd = "llhd"
        if return_nll:
            sllhd = "nll"
        muhat, sigma_mu = float("inf"), float("inf")
        dsig = self._likelihoodOfMu ( 1., return_nll=return_nll,
                allowNegativeSignals = allowNegativeSignals, corr = corr )
        muhat, sigma_mu =  dsig["muhat"], dsig["sigma_mu"]
        # llhd evaluated at mu_hat 
        if expected:
            muhat = 0.
        lmax = self.likelihood ( muhat, return_nll=return_nll )

        ret = { "muhat": muhat, "sigma_mu": sigma_mu, "lmax": lmax }
        return ret

    def _likelihoodOfMu ( self, mu : Union[float,None], 
            return_nll : Optional[bool] = False,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6, 
            expected : Union[bool,Text] = False ) -> float:
        """ return the likelihood, as a function of nsig
        :param mu: signal strength
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: if True, then allow muhat to become negative,
               else demand that muhat >= 0. In the presence of underfluctuations
               in the data, setting this to True results in more realistic
               approximate likelihoods.
        :param corr: correction factor:
                     ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
                     When comparing with likelihoods constructed from efficiency maps,
                     a factor of corr = 0.6 has been found to result in the best approximations.
        :returns: likelihood (float), muhat, and sigma_mu
        """
        sllhd = "llhd"
        if return_nll:
            sllhd = "nll"

        if self.upperLimitOnMu < self.expectedUpperLimitOnMu:
            ## underfluctuation. muhat = 0.
            if allowNegativeSignals:
                xa = -self.expectedUpperLimitOnMu
                xb = 1
                muhat = 0.
                if expected == False:
                    muhat = self._findMuhat( xa, xb )
                ret = self._computeLlhd(mu, muhat, return_nll = return_nll )
                return { sllhd: ret, "muhat": muhat, "sigma_mu": self.sigma_mu }
            else:
                ret = self._computeLlhd(mu, 0.0, return_nll = return_nll )
                return { sllhd: ret, "muhat": 0.0, "sigma_mu": self.sigma_mu }

        muhat = 0.
        if expected == False:
            xa = -self.expectedUpperLimitOnMu
            xb = self.expectedUpperLimitOnMu
            muhat = self._findMuhat(xa,xb)
        ret = self._computeLlhd(mu, muhat, return_nll = False )
        if return_nll:
            ret = -np.log(ret)
        return { sllhd: ret, "muhat": muhat, "sigma_mu": self.sigma_mu }

    def _getSigmaMu( self ):
        """ get the standard deviation sigma on the signal strength mu, given
        an upper limit and a central value. assumes a truncated Gaussian likelihood
        """
        # the expected scale, eq 3.24 in arXiv:1202.3415
        sigma_mu = self.expectedUpperLimitOnMu / 1.96
        #if self.newCorrectionType: # we could correct here
        #    sigma_mu = sigma_mu / ( 1. - self.corr/2.)
        return sigma_mu

    def _root_func( self, mu : float ):
        """ the root of this one should determine muhat """
        numerator = erf((self.upperLimitOnMu - mu) / self.denominator) + \
                    erf( mu / self.denominator)
        denominator = 1.0 + erf(mu / self.denominator)
        ret = numerator / denominator - self.cl
        return ret

    def _findMuhat( self, xa : float = 0., 
                        xb : Union[float,None] = None ):
        """ find muhat, in [xa,xb] 
        :param xa: lower limit of initial bracket
        :param xb: upper limit of initial bracket. if none, then max(ul,eul)
        """
        if xa == None:
            xa = 0.
        if xb == None:
            xb = max(self.upperLimitOnMu, self.expectedUpperLimitOnMu)
        c = 0
        while self._root_func(xa) * self._root_func(xb) > 0:
            xa = 2 * xa
            c += 1
            if c > 10:
                logger.error(
                    f"cannot find bracket for brent bracketing ul={self.upperLimitOnMu:.2f}, eul={self.expectedUpperLimitOnMu:.2f},r({xa:.2f})={self._root_func(xa):.2f}, r({xb:.2f})={self._root_func(xb):.2f}"
                )

        muhat = optimize.toms748(self._root_func, xa, xb, rtol=1e-07, xtol=1e-07)
        return muhat

    def _computeLlhd( self, mu, muhat, return_nll):
        if mu is None:
            mu = muhat
        sigma_mu = self.sigma_mu
        ## we could also correct here
        if self.newCorrectionType:
            sigma_mu = self.sigma_mu / ( 1. - self.corr / 2. )

        # need to account for the truncation!
        # first compute how many sigmas left of center is 0.
        Zprime = muhat / sigma_mu
        # now compute the area of the truncated gaussian
        A = stats.norm.cdf(Zprime)
        if return_nll:
            return np.log(A) - stats.norm.logpdf(mu, muhat, sigma_mu)
        return float(stats.norm.pdf(mu, muhat, sigma_mu) / A)
