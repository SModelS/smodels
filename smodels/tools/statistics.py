#!/usr/bin/env python3

"""
.. module:: statistics
   :synopsis: a module meant to collect various statistical algorithms. It contains
   the procedure that computes an approximate Gaussian likelihood from an
   expected an observer upper limit. See https://arxiv.org/abs/1202.3415.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jack Y. Araz <jack.araz@durham.ac.uk>
"""

from scipy import stats, optimize
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from scipy.special import erf
import numpy as np
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from typing import Text, Optional, Union

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

def CLsfromNLL(
    nllA: float, nll0A: float, nll: float, nll0: float, return_type: Text = "CLs-alpha"
) -> float:
    """
    compute the CLs - alpha from the NLLs
    TODO: following needs explanation
    :param nllA:
    :param nll0A:
    :param nll:
    :param nll0:
    :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs"
                        CLs-alpha: returns CLs - 0.05
                        1-CLs: returns 1-CLs value
                        CLs: returns CLs value
    :return:
    """
    assert return_type in ["CLs-alpha", "1-CLs", "CLs"], f"Unknown return type: {return_type}."
    qmu = 0.0 if 2 * (nll - nll0) < 0.0 else 2 * (nll - nll0)
    sqmu = np.sqrt(qmu)
    qA = 2 * (nllA - nll0A)
    if qA < 0.0:
        qA = 0.0
    sqA = np.sqrt(qA)
    if qA >= qmu:
        CLsb = 1.0 - stats.multivariate_normal.cdf(sqmu)
        CLb = stats.multivariate_normal.cdf(sqA - sqmu)
    else:
        CLsb = 1.0 if qA == 0.0 else 1.0 - stats.multivariate_normal.cdf((qmu + qA) / (2 * sqA))
        CLb = 1.0 if qA == 0.0 else 1.0 - stats.multivariate_normal.cdf((qmu - qA) / (2 * sqA))

    CLs = CLsb / CLb if CLb > 0 else 0.0

    if return_type == "1-CLs":
        return 1.0 - CLs
    elif return_type == "CLs":
        return CLs

    return CLs - 0.05


def determineBrentBracket(mu_hat, sigma_mu, rootfinder):
    """find a, b for brent bracketing
    :param mu_hat: mu that maximizes likelihood
    :param sigm_mu: error on mu_hat (not too reliable)
    :param rootfinder: function that finds the root (usually root_func)
    """
    sigma_mu = max(sigma_mu, 0.5)  # there is a minimum on sigma_mu
    # the root should be roughly at mu_hat + 2*sigma_mu
    a = mu_hat + 1.5 * sigma_mu
    ntrials = 20
    i = 0
    foundExtra = False
    while rootfinder(a) < 0.0:
        # if this is negative, we move it to the left
        i += 1
        a -= (i**2.0) * sigma_mu
        if i > ntrials or a < -10000.0:
            for a in [0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 0.1, -0.1, 0.01, -0.01, .001, -.001, 100., -100., .0001, -.0001 ]:
                r = rootfinder(a)
                if r > 0.0:
                    foundExtra = True
                    break
            if not foundExtra:
                logger.error(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {rootfinder(a):.2f}."
                )
                logger.error(f"mu_hat={mu_hat:.2f}, sigma_mu={sigma_mu:.2f}")
                raise SModelSError(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {rootfinder(a):.2f}."
                )
    i = 0
    foundExtra = False
    b = mu_hat + 2.5 * sigma_mu
    while rootfinder(b) > 0.0:
        # if this is positive, we move it to the right
        i += 1
        b += (i**2.0) * sigma_mu
        closestr, closest = float("inf"), None
        if i > ntrials:
            for b in [1.0, 0.0, 3.0, -1.0, 10.0, -3.0, 0.1, -0.1, -10.0, 100.0, -100.0, 1000.0, .01, -.01, .001, -.001 ]:
                root = rootfinder(b)
                if root < 0.0:
                    foundExtra = True
                    break
                if root < closestr:
                    closestr = root
                    closest = b
            if not foundExtra:
                logger.error(f"cannot find a b that is right of the root (i.e. rootfinder(b) < 0).")
                logger.error(f"closest to zero rootfinder({closest})={closestr}")
                logger.error(f"mu_hat was at {mu_hat:.2f} sigma_mu at {sigma_mu:.2f}")
                raise SModelSError()
    return a, b

def deltaChi2FromLlhd(likelihood):
    """compute the delta chi2 value from a likelihood (convenience function)"""
    if likelihood == 0.0:
        return 1e10  # a very high number is good
    elif likelihood is None:
        return None

    return -2.0 * np.log(likelihood)


def chi2FromLmax(llhd, lmax):
    """compute the chi2 from likelihood and lmax"""
    chi2 = 0.0
    if llhd > 1e-200:
        from math import log

        chi2 = 2 * log(lmax / llhd)
    if chi2 < 0.0 and llhd < 1e-100:
        # numerical inaccuracy
        chi2 = 0.0
    return chi2

