#!/usr/bin/env python3

"""
.. module:: srCombinations
   :synopsis: a module to contain the logic around combinations of signal regions
              within a single analysis, be they SL-based or pyhf-based.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, Data, UpperLimitComputer
from smodels.tools.smodelsLogging import logger
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import numpy as np
from spey import ExpectationType
from smodels.tools.speyTools import SpeyComputer

def getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=0.2, allowNegativeSignals=False):
    """
    Get combined upper limit. If covariances are given in globalInfo then
    simplified likelihood is used, else if json files are given pyhf
    cimbination is performed.

    :param nsig: list of signal events in each signal region/dataset. The list
                    should obey the ordering in globalInfo.datasetOrder.
    :param expected: return expected, not observed value
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.

    :returns: upper limit on sigma*eff
    """
    computer = SpeyComputer ( dataset, nsig )
    ret = computer.poi_upper_limit ( expected, limit_on_xsec = True )
    logger.debug("Combined upper limit : {}".format(ret))
    return ret

def getCombinedLikelihood(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0, nll=False
):
    """
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected, not observed likelihood. if "posteriori",
                     compute expected posteriori.
    :param mu: signal strength parameter mu
    """

    def likelihood(mu):
        poi_test=float(mu) if isinstance(mu, (float, int)) else mu[0]
        if expected == 'posteriori':
            return statModel.asimov_likelihood ( poi_test = poi_test, expected=ExpectationType.apriori, return_nll = nll, par_bounds = bounds, init_pars = init, **args )
        else:
            return statModel.likelihood ( poi_test = poi_test, expected=expectedDict[expected], return_nll = nll, par_bounds = bounds, init_pars = init, **args )

    lbsm = likelihood(mu)

    return lbsm

def getCombinedStatistics(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, allowNegativeSignals=False, nll=False
):
    """compute lBSM, lmax, and LSM in a single run
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected values, not observed
    """
    computer = SpeyComputer ( dataset, nsig )
    def likelihood(mu):
        poi_test=float(mu) if isinstance(mu, (float, int)) else mu[0]

        if expected == 'posteriori':
            return computer.asimov_likelihood ( poi_test = poi_test,
                    expected = "apriori", return_nll = nll )
        else:
            return computer.likelihood ( poi_test = poi_test, expected = expected,
                    return_nll = nll )

    def max_likelihood():
        # init, bounds, args = self.getSpeyInitialisation ( dataset, allowNegativeSignals, initial_bracket = False )
        if expected == 'posteriori':
            return computer.maximize_asimov_likelihood(expected="apriori" )
        else:
            return computer.maximize_likelihood(allowNegativeSignals=allowNegativeSignals, expected=expected )

    lbsm = likelihood(1.)
    lsm = likelihood(0.)
    muhat, lmax = max_likelihood()

    if lsm > lmax:
        logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lsm
        muhat = 0.0
    if lbsm > lmax:
        logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lbsm
        muhat = 1.0

    sigma_mu = computer.sigma_mu ( poi_test=muhat,expected=expected )

    return {"lbsm": lbsm, "lmax": lmax, "lsm": lsm, "muhat": muhat, "sigma_mu": sigma_mu}
