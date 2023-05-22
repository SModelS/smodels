#!/usr/bin/env python3

"""
.. module:: basicStats
   :synopsis: a module thats collects various statistical codes, that \
   are shared among other submodules.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
"""

from scipy import stats, optimize
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from scipy.special import erf
import numpy as np
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from typing import Text, Optional, Union

__all__ = [ "CLsfromNLL", "determineBrentBracket", "chi2FromLmax" ]

def CLsfromNLL(
    nllA: float, nll0A: float, nll: float, nll0: float,
    return_type: Text = "CLs-alpha" ) -> float:
    """
    compute the CLs - alpha from the NLLs
    TODO: following needs explanation
    :param nllA:
    :param nll0A:
    :param nll:
    :param nll0:
    :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs" \
                        CLs-alpha: returns CLs - 0.05 \
                        1-CLs: returns 1-CLs value \
                        CLs: returns CLs value
    :return: Cls-type value, see above
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

def determineBrentBracket(mu_hat, sigma_mu, rootfinder,
         allowNegative = True ):
    """find a, b for brent bracketing
    :param mu_hat: mu that maximizes likelihood
    :param sigm_mu: error on mu_hat (not too reliable)
    :param rootfinder: function that finds the root (usually root_func)
    :param allowNegative: if False, then do not allow a or b to become negative
    :returns: the interval a,b
    """
    sigma_mu = max(sigma_mu, 0.5)  # there is a minimum on sigma_mu
    sigma_mu = min(sigma_mu, 100.) # there is a maximum on sigma_mu
    # the root should be roughly at mu_hat + 2*sigma_mu
    a = mu_hat + 1.5 * sigma_mu
    ra = rootfinder(a)
    ntrials = 20
    i = 0
    foundExtra = False
    while ra < 0.0:
        # if this is negative, we move it to the left
        i += 1
        a -= (i**2.0) * sigma_mu
        ra = rootfinder(a)
        if i > ntrials or a < -10000.0 or ra is None or ( a < 0 and not allowNegative ):
            avalues = [0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 0.1, -0.1, 0.01, -0.01, .001, -.001, 100., -100., .0001, -.0001 ]
            if not allowNegative:
                avalues = [0.0, 1.0, 3.0, 10.0, 0.1, 0.01, .001, 100., .0001 ]
            for a in avalues:
                ra = rootfinder(a)
                if ra is None: # if cls computation failed, try with next a value
                    continue
                if ra > 0.0:
                    foundExtra = True
                    break
            if not foundExtra:
                logger.error(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {ra:.2f}."
                )
                logger.error(f"mu_hat={mu_hat:.2f}, sigma_mu={sigma_mu:.2f}")
                raise SModelSError(
                    f"cannot find an a that is left of the root. last attempt, a={a:.2f}, root = {ra:.2f}."
                )
    i = 0
    foundExtra = False
    b = mu_hat + 2.5 * sigma_mu
    rb = rootfinder(b)
    while rb > 0.0:
        # if this is positive, we move it to the right
        i += 1
        b += (i**2.0) * sigma_mu
        rb = rootfinder(b)
        if rb is None: # if cls computation failed, try with a bigger b value
            continue
        closestr, closest = float("inf"), None
        if i > ntrials or ( b < 0 and not allowNegative ):
            bvalues = [1.0, 0.0, 3.0, -1.0, 10.0, -3.0, 0.1, -0.1, -10.0, 100.0, -100.0, 1000.0, .01, -.01, .001, -.001 ]
            if not allowNegative:
                bvalues = [1.0, 0.0, 3.0, 10.0, 0.1, 100.0, 1000.0, .01, .001 ]
            for b in bvalues:
                rb = rootfinder(b)
                if rb is None: # if cls computation failed, try with next b value
                    continue
                if rb < 0.0:
                    foundExtra = True
                    break
                if rb < closestr:
                    closestr = rb
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
    if llhd is None or lmax is None:
        return None
    chi2 = 0.0
    if llhd > 1e-200:
        from math import log

        chi2 = 2 * log(lmax / llhd)
    if chi2 < 0.0 and llhd < 1e-100:
        # numerical inaccuracy
        chi2 = 0.0
    return chi2
