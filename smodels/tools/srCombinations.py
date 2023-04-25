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
from smodels.tools.statsTools import StatsComputer
# from test.debug import printTo

def getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=0.2):
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
    computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
           marginalize = dataset._marginalize, normalize_sig = True )
    ret = computer.poi_upper_limit ( expected = expected, limit_on_xsec = True )
    logger.debug("pyhf upper limit : {}".format(ret))
    return ret

def getCombinedLikelihood(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0
):
    """compute only lBSM
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected, not observed likelihood. if "posteriori",
                     compute expected posteriori.
    :param mu: signal strength parameter mu
    """
    if dataset.type == "pyhf":
        # Getting the path to the json files
        # Loading the jsonFiles
        computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
               marginalize = dataset._marginalize, normalize_sig = False )
        lbsm = computer.likelihood ( poi_test = mu, expected = expected, return_nll = False )
        return lbsm
    lbsm = getCombinedSimplifiedLikelihood(
        dataset, nsig, marginalize, deltas_rel, expected=expected, mu=mu )
    return lbsm

def getCombinedPyhfStatistics(
    dataset, nsig, marginalize, deltas_rel, nll=False, expected=False,
    allowNegativeSignals=False
):
        # Getting the path to the json files
        # Loading the jsonFiles
        computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
               marginalize = dataset._marginalize, normalize_sig = False )
        ret = computer.get_five_values ( expected = expected, allowNegativeSignals =allowNegativeSignals, return_nll = False )
        return ret

def getCombinedStatistics(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, allowNegativeSignals=False
):
    """compute lBSM, lmax, and LSM in a single run
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected values, not observed
    """
    if dataset.type == "pyhf":
        return getCombinedPyhfStatistics ( dataset, nsig, marginalize, deltas_rel,
            deltas_rel, expected=expected, allowNegativeSignals=allowNegativeSignals)
    cslm = getCombinedSimplifiedStatistics( dataset, nsig, marginalize,
        deltas_rel, expected=expected, allowNegativeSignals=allowNegativeSignals,
    )
    return cslm

def getCombinedSimplifiedLikelihood(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0
):
    """
    Computes the combined simplified likelihood to observe nobs events, given a
    predicted signal "nsig", with nsig being a vector with one entry per
    dataset.  nsig has to obey the datasetOrder. Deltas is the error on
    the signal.
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected likelihood, not observed
    :param mu: signal strength parameter mu
    :returns: likelihood to observe nobs events (float)
    """
    # computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
    #       marginalize = dataset._marginalize, normalize_sig = True )
    # ret = computer.likelihood ( poi_test = mu, expected = expected )
    for k, v in enumerate(nsig):
        nsig[k] = v * mu

    if dataset.type != "simplified":
        logger.error(
            "Asked for combined simplified likelihood, but no covariance given: %s" % dataset.type
        )
        return None
    if len(dataset.origdatasets) == 1:
        if isinstance(nsig, list):
            nsig = nsig[0]
        return dataset.origdatasets[0].likelihood(nsig, marginalize=marginalize)
    bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
    nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
    cov = dataset.globalInfo.covariance
    computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
    computer.transform ( expected )
    ret = computer.likelihood(1., marginalize=marginalize)
    return ret

def getCombinedSimplifiedStatistics(
    dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False
):
    """compute likelihood at maximum, for simplified likelihoods only"""
    if dataset.type != "simplified":
        return {"lmax": -1.0, "muhat": None, "sigma_mu": None}
    computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
           marginalize = dataset._marginalize, normalize_sig = False )
    ret = computer.get_five_values ( expected = expected,
            allowNegativeSignals = allowNegativeSignals, check_for_maxima = True )
    return ret
