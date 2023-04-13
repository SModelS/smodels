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
from smodels.tools.speyTools import getSpeyInitialisation

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
    expectedDict = {False:ExpectationType.observed,
                    True:ExpectationType.apriori,
                    "posteriori":ExpectationType.aposteriori}
    if expected not in expectedDict.keys():
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    if dataset.type in ["simplified","pyhf"]:
        #statModel = dataset.statModel ###For future API
        statModel = dataset.getStatModel(nsig)
        config = statModel.backend.model.config()
        init, bounds, args = getSpeyInitialisation ( dataset, False )
        #options = { "maxiter": 1000, "method": "SLSQP", "ntrials": 3 }
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
        ret = mu_ul*statModel.xsection
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
        init[1:] = statModel.backend.model.observed - statModel.backend.model.background
        for i in range ( len ( statModel.backend.model.observed ) ):
            x = np.sqrt ( statModel.backend.model.covariance[i][i] )
            bounds[i+1]=(-7*x,7*x)

    #print ( "lbsm init pars in srCombinations are", init[:3] )
    #print ( "lbsm bounds    in srCombinations are", bounds[:3] )

    def likelihood(mu):
        poi_test=float(mu) if isinstance(mu, (float, int)) else mu[0]
        init, bounds, args = getSpeyInitialisation ( dataset, allowNegativeSignals)
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
    elif dataset.type == "simplified":
        if type(nsig)==tuple:
            nsig = np.array(nsig)
    else:
        return {"lmax": -1.0, "muhat": None, "sigma_mu": None}

    args={}

    statModel = dataset.getStatModel(nsig)

    # statModel = dataset.statModel ###For future API

    config = statModel.backend.model.config()
    init = config.suggested_init
    bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds] if dataset.type=='simplified' else config.suggested_bounds
    if allowNegativeSignals:
        bounds[config.poi_index] = (config.minimum_poi, 100)
    else:
        bounds[config.poi_index] = (0, 100)
    if dataset.type == "simplified":
        init[1:] = statModel.backend.model.observed - statModel.backend.model.background
        for i in range ( len ( statModel.backend.model.observed ) ):
            x = np.sqrt ( statModel.backend.model.covariance[i][i] )
            bounds[i+1]=(-7*x,7*x)

    # print ( "for muhat init pars in srCombinations are", init[:3] )
    # print ( "for muhat bounds    in srCombinations are", bounds[:3] )

    def likelihood(mu):
        poi_test=float(mu) if isinstance(mu, (float, int)) else mu[0]
        init, bounds, args = getSpeyInitialisation ( dataset, allowNegativeSignals)
        if expected == 'posteriori':
            return statModel.asimov_likelihood ( poi_test = poi_test, expected=ExpectationType.apriori, return_nll = nll, par_bounds = bounds, init_pars = init, **args )
        else:
            return statModel.likelihood ( poi_test = poi_test, expected=expectedDict[expected], return_nll = nll, par_bounds = bounds, init_pars = init, **args )

    def max_likelihood():
        init, bounds, args = getSpeyInitialisation ( dataset, allowNegativeSignals)
        if expected == 'posteriori':
            return statModel.maximize_asimov_likelihood(expected=ExpectationType.apriori, test_statistics="qmutilde", return_nll=nll, par_bounds=bounds, init_pars=init, **args )
        else:
            return statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds, init_pars = init, **args )


    lbsm = likelihood(1.)
    lsm = likelihood(0.)
    muhat, lmax = max_likelihood()

    while abs(muhat - bounds[config.poi_index][1]) <= 0.1:
        logger.debug('Muhat reached the upper bound. Will try again after increasing the upper bound.')
        bounds[config.poi_index] = (bounds[config.poi_index][1], bounds[config.poi_index][1]*10)
        muhat, lmax = max_likelihood()

    if lsm > lmax:
        logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lsm
        muhat = 0.0
    if lbsm > lmax:
        logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lbsm
        muhat = 1.0

    # sys.exit()
    test_statistics = "q" if allowNegativeSignals else "qmutilde"
    sigma_mu = statModel.sigma_mu(poi_test=muhat,expected=expectedDict[expected],test_statistics=test_statistics)

    return {"lbsm": lbsm, "lmax": lmax, "lsm": lsm, "muhat": muhat, "sigma_mu": sigma_mu}

    # if dataset.type == "pyhf":
        # return getCombinedPyhfStatistics ( dataset, nsig, marginalize, deltas_rel,
        #     deltas_rel, expected=expected, allowNegativeSignals=allowNegativeSignals)
    # cslm = getCombinedSimplifiedStatistics( dataset, nsig, marginalize,
    #     deltas_rel, expected=expected, allowNegativeSignals=allowNegativeSignals,
    # )
    # return cslm

def _getPyhfComputer(dataset, nsig, normalize=True):
    """create the pyhf ul computer object
    :param normalize: if true, normalize nsig
    :returns: pyhf upper limit computer, and combinations of signal regions
    """
    # Getting the path to the json files
    jsonFiles = [js for js in dataset.globalInfo.jsonFiles]
    jsons = dataset.globalInfo.jsons.copy()
    # datasets = [ds.getID() for ds in dataset._datasets]
    datasets = [ds.getID() for ds in dataset.origdatasets]
    total = sum(nsig)
    if total == 0.0:  # all signals zero? can divide by anything!
        total = 1.0
    if normalize:
        nsig = [
            s / total for s in nsig
        ]  # Normalising signals to get an upper limit on the events count
    # Filtering the json files by looking at the available datasets
    for jsName in dataset.globalInfo.jsonFiles:
        if all([ds not in dataset.globalInfo.jsonFiles[jsName] for ds in datasets]):
            # No datasets found for this json combination
            jsIndex = jsonFiles.index(jsName)
            jsonFiles.pop(jsIndex)
            jsons.pop(jsIndex)
            continue
        if not all([ds in datasets for ds in dataset.globalInfo.jsonFiles[jsName]]):
            # Some SRs are missing for this json combination
            logger.error( "Wrong json definition in globalInfo.jsonFiles for json : %s" % jsName)
    logger.debug("list of datasets: {}".format(datasets))
    logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
    # Constructing the list of signals with subsignals matching each json
    nsignals = list()
    for jsName in jsonFiles:
        subSig = list()
        for srName in dataset.globalInfo.jsonFiles[jsName]:
            try:
                index = datasets.index(srName)
            except ValueError:
                line = (
                    f"{srName} signal region provided in globalInfo is not in the list of datasets, {jsName}:{','.join(datasets)}"
                )
                raise ValueError(line)
            sig = nsig[index]
            subSig.append(sig)
        nsignals.append(subSig)
    # Loading the jsonFiles, unless we already have them (because we pickled)
    from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer

    data = PyhfData(nsignals, jsons, jsonFiles)
    if data.errorFlag:
        return None
    if hasattr(dataset.globalInfo, "includeCRs"):
        includeCRs = dataset.globalInfo.includeCRs
    else:
        includeCRs = False
    ulcomputer = PyhfUpperLimitComputer(data, includeCRs=includeCRs,
                                        lumi=dataset.getLumi() )
    # _pyhfcomputers[idt] = ulcomputer
    return ulcomputer

# !TP - not used anymore, merged in to getCombinedLikelihood()
# def getCombinedSimplifiedLikelihood(dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0):
#     """
#     Computes the combined simplified likelihood to observe nobs events, given a
#     predicted signal "nsig", with nsig being a vector with one entry per
#     dataset.  nsig has to obey the datasetOrder. Deltas is the error on
#     the signal.
#     :param nsig: predicted signal (list)
#     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
#     :param expected: compute expected likelihood, not observed
#     :param mu: signal strength parameter mu
#     :returns: likelihood to observe nobs events (float)
#     """
#     if dataset.type != "simplified":
#         logger.error(
#             "Asked for combined simplified likelihood, but no covariance given: %s" % dataset.type
#         )
#         return None
#
#     args={"marginalize":marginalize}
#
#     statModel = dataset.getStatModel(nsig)
#
#     expectedDict = {False:ExpectationType.observed,
#                     True:ExpectationType.apriori,
#                     "posteriori":ExpectationType.aposteriori}
#     if expected not in expectedDict.keys():
#         logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
#         return None
#
#     config = statModel.backend.model.config()
#     bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
#     if allowNegativeSignals:
#         bounds[config.poi_index] = (config.minimum_poi, 100)
#     else:
#         bounds[config.poi_index] = (0, 100)
#
#     ret = statModel.likelihood(poi_test=mu, expected=expectedDict[expected], return_nll=False, par_bounds=bounds, **args)
#
#     return ret


# def getCombinedSimplifiedLikelihood(dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0):
#     """
#     Computes the combined simplified likelihood to observe nobs events, given a
#     predicted signal "nsig", with nsig being a vector with one entry per
#     dataset.  nsig has to obey the datasetOrder. Deltas is the error on
#     the signal.
#     :param nsig: predicted signal (list)
#     :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
#     :param expected: compute expected likelihood, not observed
#     :param mu: signal strength parameter mu
#     :returns: likelihood to observe nobs events (float)
#     """
#     for k, v in enumerate(nsig):
#         nsig[k] = v * mu
#
#     if dataset.type != "simplified":
#         logger.error(
#             "Asked for combined simplified likelihood, but no covariance given: %s" % dataset.type
#         )
#         return None
#     if len(dataset.origdatasets) == 1:
#         if isinstance(nsig, list):
#             nsig = nsig[0]
#         return dataset.origdatasets[0].likelihood(nsig, marginalize=marginalize)
#     bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
#     nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
#     if expected == True:
#         nobs = [x.dataInfo.expectedBG for x in dataset.origdatasets]
#     if expected == False:
#         nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
#     cov = dataset.globalInfo.covariance
#     computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
#     if expected == "posteriori":
#         theta_hat, _ = computer.findThetaHat ( 0. )
#         nobs = [float(x + y) for x, y in zip(bg, theta_hat)]
#         computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
#     return computer.likelihood(1., marginalize=marginalize)

#!TP not used anymore, merged into getCombinedStatistics()
# def getCombinedSimplifiedStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
#     """compute likelihood at maximum, for simplified likelihoods only"""
#     if dataset.type != "simplified":
#         return {"lmax": -1.0, "muhat": None, "sigma_mu": None}
#     args={"marginalize":marginalize}
#
#     if type(nsig)==tuple:
#         nsig = np.array(nsig)
#
#     statModel = dataset.getStatModel(nsig)
#
#     expectedDict = {False:ExpectationType.observed,
#                     True:ExpectationType.apriori,
#                     "posteriori":ExpectationType.aposteriori}
#     if expected not in expectedDict.keys():
#         logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
#         return None
#
#     config = statModel.backend.model.config()
#     bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
#     if allowNegativeSignals:
#         bounds[config.poi_index] = (config.minimum_poi, 100)
#     else:
#         bounds[config.poi_index] = (0, 100)
#
#     muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds, **args)
#     while abs(muhat - bounds[config.poi_index][1]) <= 0.1:
#         logger.debug('Muhat reached the upper bound. Will try again after increasing the upper bound.')
#         bounds[config.poi_index] = (bounds[config.poi_index][1], bounds[config.poi_index][1]*10)
#         muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds, **args)
#
#     lbsm = statModel.likelihood ( poi_test = 1., expected=expectedDict[expected], return_nll = nll, **args )
#     lsm = statModel.likelihood ( poi_test = 0., expected=expectedDict[expected], return_nll = nll, **args )
#     if lsm > lmax:
#         logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         lmax = lsm
#         muhat = 0.0
#     if lbsm > lmax:
#         logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         lmax = lbsm
#         muhat = 1.0
#     test_statistics = "q" if allowNegativeSignals else "qmutilde"
#     sigma_mu = statModel.sigma_mu(poi_test=muhat,expected=expectedDict[expected],test_statistics=test_statistics)
#     return {"muhat": muhat, "sigma_mu": sigma_mu, "lmax": lmax, "lbsm": lbsm, "lsm": lsm }


# def getCombinedSimplifiedStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
#     """compute likelihood at maximum, for simplified likelihoods only"""
#     if dataset.type != "simplified":
#         return {"lmax": -1.0, "muhat": None, "sigma_mu": None}
#     nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
#     bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
#     if expected == True:
#         # nobs = [ x.dataInfo.expectedBG for x in dataset._datasets]
#         nobs = [x.dataInfo.expectedBG for x in dataset.origdatasets]
#         # nobs = [int(np.round(x.dataInfo.expectedBG)) for x in dataset._datasets]
#     bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
#     cov = dataset.globalInfo.covariance
#     if type(nsig) in [list, tuple]:
#         nsig = np.array(nsig)
#     computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
#     if expected == "posteriori":
#         theta_hat, _ = computer.findThetaHat ( 0. )
#         nobs = [float(x + y) for x, y in zip(bg, theta_hat)]
#         computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
#     ret = computer.findMuHat(allowNegativeSignals=allowNegativeSignals, extended_output=True)
#     lbsm = computer.likelihood ( 1., marginalize = marginalize )
#     lsm = computer.likelihood ( 0., marginalize = marginalize )
#     lmax = ret["lmax"]
#     if lsm > lmax:
#         muhat = ret["muhat"]
#         logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         ret["lmax"] = lsm
#         ret["muhat"] = 0.0
#     if lbsm > lmax:
#         muhat = ret["muhat"]
#         logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
#         ret["lmax"] = lbsm
#         ret["muhat"] = 1.0
#     ret["lbsm"] = lbsm
#     ret["lsm"] = lsm
#     return ret
