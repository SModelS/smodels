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
    #expectedDict = {False:ExpectationType.observed,
    #                True:ExpectationType.apriori,
    #                "posteriori":ExpectationType.aposteriori}
    #if expected not in expectedDict.keys():
    #    logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
    #    return None

    if dataset.type in ["simplified","pyhf"]:
        #statModel = dataset.statModel ###For future API
        computer = SpeyComputer ( dataset, nsig )
        mu_ul = computer.poi_upper_limit ( expected )
        """
        statModel = dataset.getStatModel(nsig)
        config = statModel.backend.model.config()
        init, bounds, args = getSpeyInitialisation ( dataset, False, False )
        try:
            mu_ul = statModel.poi_upper_limit(expected=expectedDict[expected], par_bounds=bounds, init_pars = init, **args )
        except ValueError as e: # try with SLSQP or BFGS  and different bracket
            logger.warning ( f"when computing upper limit for SL: {e}. Will try with other method" )
            alternateMethod ( args )
            #if "xrtol" in args:
            #    args.pop ( "xrtol" )
            mu_ul = statModel.poi_upper_limit(expected=expectedDict[expected], par_bounds=bounds, init_pars = init, **args )
        #mu_ul = statModel.poi_upper_limit(expected=expectedDict[expected],par_bounds=bounds, init_pars = init, **options )
        if False:
            print ( "in srCombinations mu_ul is", mu_ul )
            muhat = statModel.maximize_likelihood ( par_bounds = bounds, init_pars  = init )
            print ( "in srCombinations muhat is", muhat )
            # import sys; sys.exit()
        while abs(mu_ul - bounds[config.poi_index][1]) <= 0.1:
            logger.debug('Upper limit on poi reached the upper bound. Will try again after increasing the upper bound.')
            bounds[config.poi_index] = (bounds[config.poi_index][1], bounds[config.poi_index][1]*10)
            mu_ul = statModel.poi_upper_limit(expected=ExpectationType.apriori,par_bounds=bounds, init_pars = init )
        """
        ret = mu_ul*computer.xsection
        logger.debug("Combined upper limit : {}".format(ret))
        return ret
    else:
        logger.error(
            "no covariance matrix or json file given in globalInfo.txt for %s"
            % dataset.globalInfo.id
        )
        raise SModelSError(
            "no covariance matrix or json file given in globalInfo.txt for %s"
            % dataset.globalInfo.id
        )


# def getCombinedUpperLimitFor(dataset, nsig, expected=False, deltas_rel=0.2):
#     """
#     Get combined upper limit. If covariances are given in globalInfo then
#     simplified likelihood is used, else if json files are given pyhf
#     cimbination is performed.
#
#     :param nsig: list of signal events in each signal region/dataset. The list
#                     should obey the ordering in globalInfo.datasetOrder.
#     :param expected: return expected, not observed value
#     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
#
#     :returns: upper limit on sigma*eff
#     """
#
#     if dataset.type == "simplified":
#         cov = dataset.globalInfo.covariance
#         if type(cov) != list:
#             raise SModelSError("covariance field has wrong type: %s" % type(cov))
#         if len(cov) < 1:
#             raise SModelSError("covariance matrix has length %d." % len(cov))
#
#         computer = UpperLimitComputer(ntoys=10000)
#
#         nobs = [x.dataInfo.observedN for x in dataset._datasets]
#         bg = [x.dataInfo.expectedBG for x in dataset._datasets]
#         no = nobs
#
#         d = Data(
#             observed=no,
#             backgrounds=bg,
#             covariance=cov,
#             third_moment=None,
#             nsignal=nsig,
#             deltas_rel=deltas_rel,
#             lumi=dataset.getLumi(),
#         )
#         ret = computer.getUpperLimitOnSigmaTimesEff(d, marginalize=dataset._marginalize, expected=expected)
#         logger.debug("SL upper limit : {}".format(ret))
#         return ret
#     elif dataset.type == "pyhf":
#         logger.debug("Using pyhf")
#         if all([s == 0 for s in nsig]):
#             logger.warning("All signals are empty")
#             return None
#         ulcomputer = _getPyhfComputer(dataset, nsig)
#         ret = ulcomputer.getUpperLimitOnSigmaTimesEff(expected=expected)
#         logger.debug("pyhf upper limit : {}".format(ret))
#         return ret
#     else:
#         logger.error(
#             "no covariance matrix or json file given in globalInfo.txt for %s"
#             % dataset.globalInfo.id
#         )
#         raise SModelSError(
#             "no covariance matrix or json file given in globalInfo.txt for %s"
#             % dataset.globalInfo.id
#         )


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
    expectedDict = {False:ExpectationType.observed,
                    True:ExpectationType.apriori,
                    "posteriori":ExpectationType.aposteriori}
    if expected not in expectedDict.keys():
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    if dataset.type == "pyhf":
        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for pyhf backend.")
        if marginalize == True:
            logger.error('Pyhf backend cannot marginalize likelihood.')
            return None

    statModel = dataset.getStatModel(nsig)

    # statModel = dataset.statModel ###For future API

    config = statModel.backend.model.config()
    if mu < config.minimum_poi:
        logger.error ( f'Calling likelihood for {dataset.globalInfo.id} (using combination of SRs) for a mu giving a negative total yield. mu = {mu} and minimum_mu = {config.minimum_poi}.' )
        return None
    bounds = config.suggested_bounds
    bounds[config.poi_index] = (max(mu-0.1,config.minimum_poi),mu+0.1)
    init = config.suggested_init
    args={}
    if dataset.type=="simplified":
        assert config.poi_index == 0, f"Error: I assume the poi index to be zero, not {config.poi_index}"

    #print ( "lbsm init pars in srCombinations are", init[:3] )
    #print ( "lbsm bounds    in srCombinations are", bounds[:3] )

    def likelihood(mu):
        poi_test=float(mu) if isinstance(mu, (float, int)) else mu[0]
        init, bounds, args = getSpeyInitialisation ( dataset, allowNegativeSignals, initial_bracket = False )
        if expected == 'posteriori':
            return statModel.asimov_likelihood ( poi_test = poi_test, expected=ExpectationType.apriori, return_nll = nll, par_bounds = bounds, init_pars = init, **args )
        else:
            return statModel.likelihood ( poi_test = poi_test, expected=expectedDict[expected], return_nll = nll, par_bounds = bounds, init_pars = init, **args )

    lbsm = likelihood(mu)

    return lbsm

# !TP not used anymore, merged into getCombinedStatistics()
# def getCombinedPyhfStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
#     expectedDict = {False:ExpectationType.observed,
#                     True:ExpectationType.apriori,
#                     "posteriori":ExpectationType.aposteriori}
#     if expected not in expectedDict.keys():
#         logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
#         return None
#
#     if deltas_rel != 0.2:
#         logger.warning("Relative uncertainty on signal not supported by spey for pyhf backend.")
#     if marginalize == True:
#         logger.error('Pyhf backend cannot marginalize likelihood.')
#
#     statModel = dataset.getStatModel(nsig)
#
#     config = statModel.backend.model.config()
#     bounds = config.suggested_bounds
#     if allowNegativeSignals:
#         bounds[config.poi_index] = (config.minimum_poi, 100)
#     else:
#         bounds[config.poi_index] = (0, 100)
#
#     muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds)
#     while abs(muhat - bounds[config.poi_index][1]) <= 0.1:
#         logger.debug('Muhat reached the upper bound. Will try again after increasing the upper bound.')
#         bounds[config.poi_index] = (bounds[config.poi_index][1], bounds[config.poi_index][1]*10)
#         muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds)
#
#     lbsm = statModel.likelihood ( poi_test = 1., expected=expectedDict[expected], return_nll = nll)
#     lsm = statModel.likelihood ( poi_test = 0., expected=expectedDict[expected], return_nll = nll)
#     if lsm > lmax:
#         logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         lmax = lsm
#         muhat = 0.0
#     if lbsm > lmax:
#         logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         lmax = lbsm
#         muhat = 1.0
#
#     test_statistics = "q" if allowNegativeSignals else "qmutilde"
#     sigma_mu = statModel.sigma_mu(poi_test=muhat,expected=expectedDict[expected],test_statistics=test_statistics)
#
#     return {"lbsm": lbsm, "lmax": lmax, "lsm": lsm, "muhat": muhat, "sigma_mu": sigma_mu}


# def getCombinedPyhfStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
#         # Getting the path to the json files
#         # Loading the jsonFiles
#         ulcomputer = _getPyhfComputer(dataset, nsig, False)
#         index = ulcomputer.getBestCombinationIndex()
#         lbsm = ulcomputer.likelihood(mu=1.0, workspace_index=index, expected=expected)
#         lmax = ulcomputer.lmax(
#             workspace_index=index, expected=expected, allowNegativeSignals=allowNegativeSignals
#         )
#         muhat = None
#         try:
#             muhat = float(ulcomputer.muhat)
#         except AttributeError:
#             pass
#         sigma_mu = ulcomputer.sigma_mu
#         ulcomputer = _getPyhfComputer(dataset, [0.0] * len(nsig), False)
#         lsm = ulcomputer.likelihood(mu=0.0, workspace_index=index, expected=expected)
#         return {"lbsm": lbsm, "lmax": lmax, "lsm": lsm, "muhat": muhat, "sigma_mu": sigma_mu}

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
