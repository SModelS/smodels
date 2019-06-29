#!/usr/bin/env python3

"""
.. module:: statistics
   :synopsis: a module meant to collect various statistical algorithms. 
              For now it only contains the procedure that computes an
              approximate Gaussian likelihood from an expected an observer upper
              limit

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from scipy import stats, optimize, integrate, special
from scipy.special import erf
import numpy as np

def likelihoodFromLimits( upperLimit, expectedUpperLimit, nsig ):
    """ computes the likelihood from an expected and an observed upper limit.
    :param upperLimit: observed upper limit, as a yield (i.e. unitless)
    :param expectedUpperLimit: expected upper limit, also as a yield
    :param nSig: number of signal events

    :returns: likelihood
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
    likelihood = stats.norm.pdf ( nsig, mumax, sigma_exp )
    # print ( ">>> mu=%.2f, mumax=%.3f llhd = %.5g" % ( mu, mumax, likelihood ),"ul=",ul,"eul=",eul )
    return likelihood

def chi2FromLlhd ( likelihood ):
    """ compute the chi2 value from a likelihood (convenience function) """
    return -2. * np.log ( likelihood )

