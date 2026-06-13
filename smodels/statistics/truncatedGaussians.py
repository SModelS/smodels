#!/usr/bin/env python3

"""
.. module:: truncatedGaussian
   :synopsis: a module that contains the code that computes an approximate
              Gaussian likelihood from an evaluationType an observer upper limit. See
              https://arxiv.org/abs/1202.3415.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "TruncatedGaussians" ]

from scipy import stats, optimize
from smodels.base.smodelsLogging import logger
from scipy.special import erf
import numpy as np
from smodels.statistics.basicStats import observed, NllEvalType
from typing import  Optional, Union, Dict

class TruncatedGaussians:
    """ likelihood computer based on the trunacated Gaussian approximation, see
         arXiv:1202.3415 """
    __slots__ = [ "upperLimitOnMu", "expectedUpperLimitOnMu", "corr",
                  "sigma_mu", "denominator", "cl", "allowNegativeSignals" ]

    # the correction we introduced a long time ago was have exclusions
    # in mind only. for discovery mode we need to correct differently
    newCorrectionType = False

    def __init__  ( self, upperLimitOnMu : float, expectedUpperLimitOnMu : float,
                    corr : Optional[float] = 0.6, cl=.95 ):
        """
        :param upperLimitOnMu: observed upper limit on signal strength mu
        :param expectedUpperLimitOnMu: evaluationType upper limit
        on signal strength mu
        :param corr: correction factor:
        ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
        When comparing with likelihoods constructed from efficiency maps,
        a factor of corr = 0.6 has been found to result
        in the best approximations.
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
        self.sigma_mu = self._getSigmaMu()  # the evaluationType scale, eq 3.24 in arXiv:1202.3415
        self.denominator = np.sqrt(2.0) * self.sigma_mu
        self.allowNegativeSignals = False
        self.cl = cl

    """
    def likelihood ( self, mu : Union[float,None], return_nll : Optional[bool]=False,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            evaluationType : NllEvalType = observed ) -> Union[None,float]:
        ret = self.nll( mu, allowNegativeSignals,
                corr, evaluationType )
        if return_nll:
            return ret
        return math.exp ( - ret )
    """

    def nll( self, mu : Union[float,None],
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            evaluationType : NllEvalType = observed ) -> Union[None,float]:
        """ return the nll, as a function of mu

        :param mu: number of signal events, if None then mu = muhat
        :param allowNegativeSignals: if True, then allow muhat to become negative,\
               else demand that muhat >= 0. In the presence of underfluctuations\
               in the data, setting this to True results in more realistic\
               approximate likelihoods.

        :returns: nll (float)
        """
        dsig = self._nllOfMu ( mu, allowNegativeSignals = allowNegativeSignals,
                corr = corr, evaluationType = evaluationType )
        ret = dsig["nll"]
        return ret

    def nll_min( self, return_nll : Optional[bool]=False,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            evaluationType : NllEvalType = observed ) -> Dict:
        """
        Return the nll at muhat

        :param mu: number of signal events, if None then mu = muhat
        :param allowNegativeSignals: if True, then allow muhat to become negative,
        else demand that muhat >= 0. In the presence of underfluctuations
        in the data, setting this to True results in more realistic
        approximate likelihoods.

        :returns: dictionary with likelihood (float), muhat, and sigma_mu
        """

        muhat, sigma_mu = float("inf"), float("inf")
        dsig = self._nllOfMu ( 1.,
                allowNegativeSignals = allowNegativeSignals, corr = corr )
        muhat, sigma_mu =  dsig["muhat"], dsig["sigma_mu"]
        # llhd evaluated at mu_hat
        if evaluationType != observed:
            muhat = 0.
        nll_min = self.nll ( muhat )

        ret = { "muhat": muhat, "sigma_mu": sigma_mu, "nll_min": nll_min }
        return ret

    def _nllOfMu ( self, mu : float,
            allowNegativeSignals : Optional[bool] = True,
            corr : Optional[float] = 0.6,
            evaluationType : NllEvalType = observed ) -> Dict:
        """ return the nll, as a function of nsig

        :param mu: signal strength
        :param allowNegativeSignals: if True, then allow muhat to become negative,
        else demand that muhat >= 0. In the presence of underfluctuations
        in the data, setting this to True results in more realistic
        approximate likelihoods.
        :param corr: correction factor:
        ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
        When comparing with likelihoods constructed from efficiency maps,
        a factor of corr = 0.6 has been found to result in the best approximations.
        :returns: nll (float), muhat, and sigma_mu
        """
        if self.upperLimitOnMu < self.expectedUpperLimitOnMu:
            ## underfluctuation. muhat = 0.
            if allowNegativeSignals:
                xa = -self.expectedUpperLimitOnMu
                xb = 1
                muhat = 0.
                if evaluationType == observed:
                    muhat = self._findMuhat( xa, xb )
                ret = self._computeNLL(mu, muhat )
                return { "nll": ret, "muhat": muhat, "sigma_mu": self.sigma_mu }
            else:
                ret = self._computeNLL(mu, 0.0 )
                return { "nll": ret, "muhat": 0.0, "sigma_mu": self.sigma_mu }

        muhat = 0.
        if evaluationType == observed:
            xa = -self.expectedUpperLimitOnMu
            xb = self.expectedUpperLimitOnMu
            muhat = self._findMuhat(xa,xb)
        ret = self._computeNLL(mu, muhat )
        return { "nll": ret, "muhat": muhat, "sigma_mu": self.sigma_mu }

    def _getSigmaMu( self ) -> float :
        """ get the standard deviation sigma on the signal strength mu, given
        an upper limit and a central value.
        assumes a truncated Gaussian likelihood
        """
        # the evaluationType scale, eq 3.24 in arXiv:1202.3415
        sigma_mu = self.expectedUpperLimitOnMu / 1.96
        #if self.newCorrectionType: # we could correct here
        #    sigma_mu = sigma_mu / ( 1. - self.corr/2.)
        return sigma_mu

    def _root_func( self, mu : float ) -> float:
        """ the root of this one should determine muhat """
        numerator = erf((self.upperLimitOnMu - mu) / self.denominator) + \
                    erf( mu / self.denominator)
        denominator = 1.0 + erf(mu / self.denominator)
        if denominator == 0. and np.isclose ( numerator, denominator ):
            ## 0 / 0 = 1 (in this case)
            return 1. - self.cl
        ret = numerator / denominator - self.cl
        return ret

    def _findMuhat( self, xa : float = 0.,
                        xb : Optional[float] = None ) -> float:
        """ find muhat, in [xa,xb]

        :param xa: lower limit of initial bracket
        :param xb: upper limit of initial bracket. if none, then max(ul,eul)
        """
        if xa == None:
            xa = 0.
        if xb == None:
            xb = max(self.upperLimitOnMu, self.expectedUpperLimitOnMu)
        c = 0
        ra, rb = self._root_func(xa), self._root_func(xb)
        while ra * rb > 0:
            c += 1
            if c > 10:
                line = f"cannot find bracket for brent bracketing ul={self.upperLimitOnMu:.2f}, eul={self.expectedUpperLimitOnMu:.2f},r({xa:.2f})={self._root_func(xa):.2f}, r({xb:.2f})={self._root_func(xb):.2f}"
                logger.error( line )
                raise ValueError ( line )
            while ra < 0.:
                xa = xa - 1.5 * abs(xa)-.3
                ra = self._root_func ( xa )
            while rb > 0.:
                xb = xb + 1.5 * abs(xb)+.3
                rb = self._root_func ( xb )

        try:
            muhat,rootResults = optimize.toms748( self._root_func, xa, xb, rtol=1e-07,
                                                  xtol=1e-07, full_output = True  )
        except ValueError as e:
            logger.error ( f"truncated gaussian got ValueError {e}: rf({xa:.2f})={self._root_func(xa):.2f}, rf({xb:.2f})={self._root_func(xb):.2f}" )
            logger.error ( f"ul={self.upperLimitOnMu:.2f}, eul={self.expectedUpperLimitOnMu:.2f}" )
            muhat,rootResults = optimize.toms748( self._root_func, xa, xb, rtol=1e-02,
                                      xtol=1e-02, full_output = True )
            if rootResults.converged:
                logger.error ( "retry with tol=1e-2 seemed to have worked" )
            else:
                muhat,*_ = optimize.brentq(self._root_func, xa, xb, rtol=1e-02, xtol=1e-02 )
                logger.error ( "retry with brentq seemed to have worked" )
        return muhat

    def _computeNLL( self, mu : float, muhat : float ) -> float:
        if mu is None:
            mu = muhat
        sigma_mu = self.sigma_mu
        ## we could also correct here
        if self.newCorrectionType and self.corr is not None:
            sigma_mu = self.sigma_mu / ( 1. - self.corr / 2. )

        # need to account for the truncation!
        # first compute how many sigmas left of center is 0.
        Zprime = muhat / sigma_mu
        # now compute the area of the truncated gaussian
        A = stats.norm.cdf(Zprime)
        # if return_nll:
        return float ( np.log(A) - stats.norm.logpdf(mu, muhat, sigma_mu) )
        #return float(stats.norm.pdf(mu, muhat, sigma_mu) / A)

    def getUpperLimitOnMu(self, evaluationType : NllEvalType=observed, **kwargs) -> float:
        """
        Get the upper limit on the signal strength mu directly from its definition.
        """
        
        if evaluationType == observed:
            return self.upperLimitOnMu
        else:
            return self.expectedUpperLimitOnMu

    def CLs(self, **kwargs) -> None:
        """
        Dummy function in case it is needed.
        """

        return None