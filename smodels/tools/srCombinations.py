#!/usr/bin/env python3

"""
.. module:: srCombinations
   :synopsis: a module to contain the logic around combinations of signal regions
              within a single analysis, be they SL-based or pyhf-based.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

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
    computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
           marginalize = dataset._marginalize, normalize_sig = False )
    ret = computer.likelihood ( poi_test = mu, expected = expected,
                                return_nll = False )
    return ret

def getCombinedStatistics(
    dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, allowNegativeSignals=False
):
    """compute lBSM, lmax, and LSM in a single run
    :param nsig: predicted signal (list)
    :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
    :param expected: compute expected values, not observed
    """
    computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
           marginalize = dataset._marginalize, normalize_sig = False )
    ret = computer.get_five_values ( expected = expected,
            allowNegativeSignals = allowNegativeSignals, check_for_maxima = True )
    return ret

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
    computer = StatsComputer ( dataset, nsig, deltas_rel = deltas_rel,
           marginalize = dataset._marginalize, normalize_sig = False )
    ret = computer.likelihood ( poi_test = mu, expected = expected,
                                return_nll = False )
    return ret
