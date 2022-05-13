#!/usr/bin/env python3

"""
.. module:: statistics
   :synopsis: a module meant to collect various statistical algorithms. It contains
   the procedure that computes an approximate Gaussian likelihood from an
   expected an observer upper limit. See https://arxiv.org/abs/1202.3415.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from scipy import stats, optimize
from smodels.tools.smodelsLogging import logger
from scipy.special import erf
import numpy as np
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from typing import Text


def likelihoodFromLimits(
    upperLimit, expectedUpperLimit, nsig, nll=False, allowNegativeMuhat=True, corr=0.6
):
    """computes the likelihood from an expected and an observed upper limit.
    :param upperLimit: observed upper limit, as a yield (i.e. unitless)
    :param expectedUpperLimit: expected upper limit, also as a yield
    :param nsig: number of signal events, if None then nsig = mumax
    :param nll: if True, return negative log likelihood
    :param allowNegativeMuhat: if True, then allow muhat to become negative,
           else demand that muhat >= 0. In the presence of underfluctuations
           in the data, setting this to True results in more realistic
           approximate likelihoods.
    :param corr: correction factor:
                 ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
                 When comparing with likelihoods constructed from efficiency maps,
                 a factor of corr = 0.6 has been found to result in the best approximations.

    :returns: likelihood (float), muhat, and sigma_mu
    """
    assert upperLimit > 0.0

    if corr > 0.0 and upperLimit > expectedUpperLimit:
        expectedUpperLimit = expectedUpperLimit / (
            1.0 - corr * ((upperLimit - expectedUpperLimit) / (upperLimit + expectedUpperLimit))
        )

    def getSigma(ul, muhat=0.0):
        """get the standard deviation sigma, given
        an upper limit and a central value. assumes a truncated Gaussian likelihood"""
        # the expected scale, eq 3.24 in arXiv:1202.3415
        return (ul - muhat) / 1.96

    sigma_exp = getSigma(expectedUpperLimit)  # the expected scale, eq 3.24 in arXiv:1202.3415
    denominator = np.sqrt(2.0) * sigma_exp

    def root_func(x):  # we want the root of this one
        return (erf((upperLimit - x) / denominator) + erf(x / denominator)) / (
            1.0 + erf(x / denominator)
        ) - 0.95

    def find_neg_mumax(upperLimit, expectedUpperLimit, xa, xb):
        c = 0
        while root_func(xa) * root_func(xb) > 0:
            xa = 2 * xa
            c += 1
            if c > 10:
                logger.error(
                    f"cannot find bracket for brent bracketing ul={upperLimit}, eul={expectedUpperLimit},xa={xa}, xb={xb}"
                )

        mumax = optimize.brentq(root_func, xa, xb, rtol=1e-03, xtol=1e-06)
        return mumax

    def llhd(nsig, mumax, sigma_exp, nll):
        if nsig is None:
            nsig = mumax
        # need to account for the truncation!
        # first compute how many sigmas left of center is 0.
        Zprime = mumax / sigma_exp
        # now compute the area of the truncated gaussian
        A = stats.norm.cdf(Zprime)
        if nll:
            return np.log(A) - stats.norm.logpdf(nsig, mumax, sigma_exp)
        return float(stats.norm.pdf(nsig, mumax, sigma_exp) / A)

    if upperLimit < expectedUpperLimit:
        ## underfluctuation. mumax = 0.
        if allowNegativeMuhat:
            xa = -expectedUpperLimit
            xb = 1
            mumax = find_neg_mumax(upperLimit, expectedUpperLimit, xa, xb)
            return llhd(nsig, mumax, sigma_exp, nll), mumax, sigma_exp
        else:
            return llhd(nsig, 0.0, sigma_exp, nll), 0.0, sigma_exp

    fA = root_func(0.0)
    fB = root_func(max(upperLimit, expectedUpperLimit))
    if np.sign(fA * fB) > 0.0:
        ## the have the same sign
        logger.error("when computing likelihood: fA and fB have same sign")
        return None, None, None
    mumax = optimize.brentq(
        root_func, 0.0, max(upperLimit, expectedUpperLimit), rtol=1e-03, xtol=1e-06
    )
    llhdexp = llhd(nsig, mumax, sigma_exp, nll)
    return llhdexp, mumax, sigma_exp


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
    :param return_type: (Text) "CLs-alpha" or "1-CLs"
    :return:
    """
    assert return_type in ["CLs-alpha", "1-CLs"], f"Unknown return type: {return_type}."
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

    cl = 0.95
    root = CLs - 1.0 + cl
    return root


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
            for a in [0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 0.1, -0.1, 0.01, -0.01]:
                if rootfinder(a) > 0.0:
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
            for b in [1.0, 0.0, 3.0, -1.0, 10.0, -3.0, 0.1, -0.1, -10.0, 100.0, -100.0, 1000.0]:
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


def rvsFromLimits(upperLimit, expectedUpperLimit, n=1, corr=0.0):
    """
    Generates a sample of random variates, given expected and observed likelihoods.
    The likelihood is modelled as a truncated Gaussian.

    :param upperLimit: observed upper limit, as a yield (i.e. unitless)
    :param expectedUpperLimit: expected upper limit, also as a yield
    :param n: sample size
    :param corr: correction term

    :returns: sample of random variates
    """

    if corr > 0.0 and upperLimit > expectedUpperLimit:
        expectedUpperLimit = expectedUpperLimit / (
            1.0 - corr * ((upperLimit - expectedUpperLimit) / (upperLimit + expectedUpperLimit))
        )

    sigma_exp = expectedUpperLimit / 1.96  # the expected scale
    denominator = np.sqrt(2.0) * sigma_exp

    def root_func(x):  # we want the root of this one
        return (erf((upperLimit - x) / denominator) + erf(x / denominator)) / (
            1.0 + erf(x / denominator)
        ) - 0.95

    fA, fB = root_func(0.0), root_func(max(upperLimit, expectedUpperLimit))
    if np.sign(fA * fB) > 0.0:
        # the have the same sign
        logger.error("when computing likelihood: fA and fB have same sign")
        return None
    mumax = optimize.brentq(
        root_func, 0.0, max(upperLimit, expectedUpperLimit), rtol=1e-03, xtol=1e-06
    )
    ret = []
    while len(ret) < n:
        tmp = stats.norm.rvs(mumax, sigma_exp)
        if tmp > 0.0:
            ret.append(tmp)
    return ret


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


def chi2FromLimits(likelihood, upperLimit, expectedUpperLimit, corr=0.0):
    """compute the chi2 value from a likelihood (convenience function)."""
    if corr > 0.0 and upperLimit > expectedUpperLimit:
        expectedUpperLimit = expectedUpperLimit / (
            1.0 - corr * ((upperLimit - expectedUpperLimit) / (upperLimit + expectedUpperLimit))
        )
    sigma_exp = expectedUpperLimit / 1.96  # the expected scale
    l0 = 2.0 * stats.norm.logpdf(0.0, 0.0, sigma_exp)
    l = deltaChi2FromLlhd(likelihood)
    if l is None:
        return None

    return l + l0
