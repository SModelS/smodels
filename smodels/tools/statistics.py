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
from smodels.tools import runtime

def likelihoodFromLimits( upperLimit, expectedUpperLimit, nsig, nll=False ):
    """ computes the likelihood from an expected and an observed upper limit.
    :param upperLimit: observed upper limit, as a yield (i.e. unitless)
    :param expectedUpperLimit: expected upper limit, also as a yield
    :param nsig: number of signal events, if None then nsig = mumax
    :param nll: if True, return negative log likelihood

    :returns: likelihood (float)
    """
    assert ( upperLimit > 0. )

    def getSigma ( ul, muhat = 0. ):
        """ get the standard deviation sigma, given
            an upper limit and a central value. assumes a truncated Gaussian likelihood """
        # the expected scale, eq 3.24 in arXiv:1202.3415
        return ( ul - muhat ) / 1.96

    def llhd ( nsig, mumax, sigma_exp, nll ):
        if nsig == None:
            nsig = mumax
        ## need to account for the truncation!
        ## first compute how many sigmas left of center is 0.
        Zprime = mumax / sigma_exp
        ## now compute the area of the truncated gaussian
        A = stats.norm.cdf(Zprime)
        if nll:
            return np.log(A ) - stats.norm.logpdf ( nsig, mumax, sigma_exp )
        return float ( stats.norm.pdf ( nsig, mumax, sigma_exp ) / A )

    dr = 2. * ( upperLimit - expectedUpperLimit ) / ( expectedUpperLimit + upperLimit )
    if dr>runtime._drmax:
        if runtime._cap_likelihoods == False:
            logger.warning("asking for likelihood from limit but difference between oUL(%.2f) and eUL(%.2f) is too large (dr=%.2f>%.2f)" % ( upperLimit, expectedUpperLimit, dr, runtime._drmax ) )
            return None
        oldUL = upperLimit
        upperLimit = expectedUpperLimit * ( 2. + runtime._drmax ) / ( 2. - runtime._drmax )
        logger.warning("asking for likelihood from limit but difference between oUL(%.2f) and eUL(%.2f) is too large (dr=%.2f>%.2f). capping to %.2f." % \
                ( oldUL, expectedUpperLimit, dr, runtime._drmax, upperLimit ) )
        ## we are asked to cap likelihoods, so we set observed UL such that dr == drmax

    # sigma_exp = expectedUpperLimit / 1.96 # the expected scale, eq 3.24 in arXiv:1202.3415
    sigma_exp = getSigma ( expectedUpperLimit ) # the expected scale, eq 3.24 in arXiv:1202.3415
    if upperLimit < expectedUpperLimit:
        ## underfluctuation. mumax = 0.
        return llhd ( nsig, 0., sigma_exp, nll )
    def root_func ( x ): ## we want the root of this one
        return (erf((upperLimit-x)/denominator)+erf(x/denominator)) / ( 1. + erf(x/denominator)) - .95

    denominator = np.sqrt(2.) * sigma_exp

    fA = root_func ( 0. )
    fB = root_func ( max(upperLimit,expectedUpperLimit) )
    if np.sign(fA*fB) > 0.:
        ## the have the same sign
        logger.error ( "when computing likelihood: fA and fB have same sign")
        return None
    mumax = optimize.brentq ( root_func, 0., max(upperLimit, expectedUpperLimit),
                              rtol=1e-03, xtol=1e-06 )
    llhdexp = llhd ( nsig, mumax, sigma_exp, nll )
    return llhdexp

def rvsFromLimits( upperLimit, expectedUpperLimit, n=1 ):
    """
    Generates a sample of random variates, given expected and observed likelihoods.
    The likelihood is modelled as a truncated Gaussian.

    :param upperLimit: observed upper limit, as a yield (i.e. unitless)
    :param expectedUpperLimit: expected upper limit, also as a yield
    :param n: sample size

    :returns: sample of random variates
    """

    sigma_exp = expectedUpperLimit / 1.96 # the expected scale
    denominator = np.sqrt(2.) * sigma_exp
    def root_func ( x ): ## we want the root of this one
        return (erf((upperLimit-x)/denominator)+erf(x/denominator)) / ( 1. + erf(x/denominator)) - .95

    fA,fB = root_func ( 0. ), root_func ( max(upperLimit,expectedUpperLimit) )
    if np.sign(fA*fB) > 0.:
        ## the have the same sign
        logger.error ( "when computing likelihood for %s: fA and fB have same sign" % self.analysisId() )
        return None
    mumax = optimize.brentq ( root_func, 0., max(upperLimit, expectedUpperLimit), rtol=1e-03, xtol=1e-06 )
    ret = []
    while len(ret)<n:
        tmp = stats.norm.rvs ( mumax, sigma_exp )
        if tmp > 0.:
            ret.append ( tmp )
    return ret

def deltaChi2FromLlhd ( likelihood ):
    """ compute the delta chi2 value from a likelihood (convenience function) """
    if likelihood == 0.:
        return 1e10 ## a very high number is good
    elif likelihood is None:
        return None

    return -2. * np.log ( likelihood )

def chi2FromLimits ( likelihood, expectedUpperLimit ):
    """ compute the chi2 value from a likelihood (convenience function).
    """
    sigma_exp = expectedUpperLimit / 1.96 # the expected scale
    l0 = 2. * stats.norm.logpdf ( 0., 0., sigma_exp )
    l = deltaChi2FromLlhd(likelihood)
    if l is None:
        return None

    return  l + l0
