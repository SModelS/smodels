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
from spey import get_uncorrelated_region_statistical_model, get_multi_region_statistical_model, ExpectationType
from smodels.tools.physicsUnits import fb, pb

def _getBestStatModel(dataset, nsig, allow_negative_signal=False, return_mu_ul_exp_min=False):
    """
    find the index of the best expected combination.
    :param dataset: the CombinedDataSet object.
    :param nsig: list of signal yields.
    :param allow_negative_signal: if True, allow for the best fit to have negative parameter of interest, if False forbids it.
    :return: the minimal apriori expected poi upper limit, and the spey StatisticalModel object that produced that result.
    """
    # Get the list of the names of the signal regions used in the json files
    listOfSRInJson=[]
    for SRnames in dataset.globalInfo.jsonFiles.values():
        listOfSRInJson += SRnames

    patches, listOfSignals = _getPatches(dataset, nsig)
    mu_ul_exp_min = np.inf
    # Find best combination of signal regions
    for index, (patch, json) in enumerate(zip(patches,dataset.globalInfo.jsons)):
        # If the expected signal is 0 for each SR in the combined set of SRs, skip
        if all([sig==0. for sig in listOfSignals[index]]):
            continue
        # The x-section is at the level of the TheoryPrediction
        # if there are multiple sets of SRs, set a xsec_UL for the whole analysis, i.e. that uses all the SRs,
        # so that the resulting R value Is for the whole analysis
        xsec = sum(nsig)/dataset.getLumi()
        # It is possible to do differently and to set a xsec_UL on each set of SRs but that is not how it done in SModelS so far
        # xsec = sum(listOfSignals[index])/dataset.getLumi()
        statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                        signal=patch,
                                                        observed=json,
                                                        xsection=xsec
                                                        )
        # If all the SRs are used in the json files and there is only one json files, there is only one statModel.
        # No need to compute mu_ul_exp if not needed.
        if all([ds.dataInfo.dataId in listOfSRInJson for ds in dataset._datasets]) and len(dataset.globalInfo.jsons) == 1 and return_mu_ul_exp_min == False:
            return statModel
        config = statModel.backend.model.config()
        bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
        if allow_negative_signal:
            bounds[config.poi_index] = (config.minimum_poi, 100)
        else:
            bounds[config.poi_index] = (0, 100)
        mu_ul_exp = statModel.poi_upper_limit(expected=ExpectationType.apriori,allow_negative_signal=allow_negative_signal,par_bounds=bounds)
        if mu_ul_exp == None:
            continue
        elif mu_ul_exp < mu_ul_exp_min:
            mu_ul_exp_min = mu_ul_exp
            bestStatModel = statModel

    # Check if a non-combined (uncorrelated) signal region is more contraining than the best combination obtained above
    # Check if a signal region is not in the list of SR names used in the json files
    for sig,ds in zip(nsig,dataset._datasets):
        ds = ds.dataInfo
        if ds.dataId not in listOfSRInJson:
            xsec = sig/dataset.getLumi()
            # Don't bother to compute eUL again (one could do it again if needed)
            mu_ul_exp = ds.expectedUpperLimit/xsec
            if mu_ul_exp < muull_exp_min:
                logger.info("Best constraining model is a single uncorrelated model.")
                mu_ul_exp_min = mu_ul_exp
                bestStatModel = get_uncorrelated_region_statistical_model(observations=float(ds.observedN),
                                                                        backgrounds=float(ds.expectedBG),
                                                                        background_uncertainty=float(ds.bgError),
                                                                        signal_yields=float(sig),
                                                                        xsection=xsec,
                                                                        analysis=dataset.globalInfo.id,
                                                                        backend="simplified_likelihoods" # simplified likelhood backend by default
                                                                        )
    if mu_ul_exp_min == np.inf:
        logger.error(f'No minimal upper limit on POI found for {dataset.globalInfo.id}')
        return None
    if return_mu_ul_exp_min:
        return bestStatModel, mu_ul_exp_min
    else:
        return bestStatModel

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

    if dataset.type == "simplified":
        cov = dataset.globalInfo.covariance
        if type(cov) != list:
            raise SModelSError("covariance field has wrong type: %s" % type(cov))
        if len(cov) < 1:
            raise SModelSError("covariance matrix has length %d." % len(cov))

        thirdMoment = None # Need to be implemented

        nobs = [x.dataInfo.observedN for x in dataset._datasets]
        bg = [x.dataInfo.expectedBG for x in dataset._datasets]
        xsec = sum(nsig)/dataset.getLumi()

        statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                        signal=nsig,
                                                        observed=nobs,
                                                        covariance=cov,
                                                        nb=bg,
                                                        third_moment=thirdMoment,
                                                        delta_sys=deltas_rel,
                                                        xsection=xsec
                                                        )
        config = statModel.backend.model.config()
        bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
        if allowNegativeSignals:
            bounds[config.poi_index] = (config.minimum_poi, 100)
        else:
            bounds[config.poi_index] = (0, 100)
        muul = statModel.poi_upper_limit(expected=expectedDict[expected],allow_negative_signal=allowNegativeSignals,par_bounds=bounds)
        ret = muul*xsec
        logger.debug("SL upper limit : {}".format(ret))
        return ret

    elif dataset.type == "pyhf":
        logger.debug("Using pyhf")
        if all([s == 0 for s in nsig]):
            logger.warning("All signals are empty")
            return None
        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for pyhf backend.")

        if expected == True:
            statModel, mu_ul_exp_min = _getBestStatModel(dataset=dataset, nsig=nsig, allow_negative_signal=allowNegativeSignals, return_mu_ul_exp_min=True)
            return mu_ul_exp_min*statModel.xsection
        else:
            statModel = _getBestStatModel(dataset=dataset, nsig=nsig, allow_negative_signal=allowNegativeSignals, return_mu_ul_exp_min=False)
            config = statModel.backend.model.config()
            bounds = config.suggested_bounds
            if allowNegativeSignals:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            mu_ul = statModel.poi_upper_limit(expected=expectedDict[expected], allow_negative_signal=allowNegativeSignals,par_bounds=bounds)
            xsec_ul = mu_ul*statModel.xsection

            logger.debug("pyhf upper limit : {}".format(xsec_ul))
            return xsec_ul
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
        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        if expected not in expectedDict.keys():
            logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
            return None

        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for pyhf backend.")
        if marginalize == True:
            logger.error('Pyhf backend cannot marginalize likelihood.')

        statModel = _getBestStatModel(dataset, nsig)

        lbsm = statModel.likelihood(poi_test = mu, expected=expectedDict[expected], return_nll=False)

        # ulcomputer = _getPyhfComputer(dataset, nsig, False)
        # index = ulcomputer.getBestCombinationIndex()
        # lbsm = ulcomputer.likelihood(mu=mu, workspace_index=index, expected=expected)
        return lbsm
    lbsm = getCombinedSimplifiedLikelihood(
        dataset, nsig, marginalize, deltas_rel, expected=expected, mu=mu )
    return lbsm

def getCombinedPyhfStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
        # Getting the path to the json files
        # Loading the jsonFiles

        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        if expected not in expectedDict.keys():
            logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
            return None

        if deltas_rel != 0.2:
            logger.warning("Relative uncertainty on signal not supported by spey for pyhf backend.")
        if marginalize == True:
            logger.error('Pyhf backend cannot marginalize likelihood.')

        statModel = _getBestStatModel(dataset, nsig, allow_negative_signal=allowNegativeSignals)

        config = statModel.backend.model.config()
        bounds = config.suggested_bounds
        if allowNegativeSignals:
            bounds[config.poi_index] = (config.minimum_poi, 100)
        else:
            bounds[config.poi_index] = (0, 100)
        muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds)
        lbsm = statModel.likelihood ( poi_test = 1., expected=expectedDict[expected], return_nll = nll)
        lsm = statModel.likelihood ( poi_test = 0., expected=expectedDict[expected], return_nll = nll)
        if lsm > lmax:
            logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
            lmax = lsm
            muhat = 0.0
        if lbsm > lmax:
            logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
            lmax = lbsm
            muhat = 1.0

        test_statistics = "q" if allowNegativeSignals else "qmutilde"
        sigma_mu = statModel.sigma_mu(poi_test=muhat,expected=expectedDict[expected],test_statistics=test_statistics)

        return {"lbsm": lbsm, "lmax": lmax, "lsm": lsm, "muhat": muhat, "sigma_mu": sigma_mu}


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

def getCombinedSimplifiedLikelihood(dataset, nsig, marginalize=False, deltas_rel=0.2, expected=False, mu=1.0):
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
    if dataset.type != "simplified":
        logger.error(
            "Asked for combined simplified likelihood, but no covariance given: %s" % dataset.type
        )
        return None

    args={"marginalize":marginalize}
    bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
    nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
    cov = dataset.globalInfo.covariance

    thirdMoment = None
    xsec = sum(nsig)/dataset.getLumi()

    statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                    signal=nsig,
                                                    observed=nobs,
                                                    covariance=cov,
                                                    nb=bg,
                                                    third_moment=thirdMoment,
                                                    delta_sys=deltas_rel,
                                                    xsection=xsec.asNumber(pb)
                                                    )
    expectedDict = {False:ExpectationType.observed,
                    True:ExpectationType.apriori,
                    "posteriori":ExpectationType.aposteriori}
    if expected not in expectedDict.keys():
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    ret = statModel.likelihood(poi_test=mu, expected=expectedDict[expected], return_nll=False, **args)

    return ret


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

def getCombinedSimplifiedStatistics(dataset, nsig, marginalize, deltas_rel, nll=False, expected=False, allowNegativeSignals=False):
    """compute likelihood at maximum, for simplified likelihoods only"""
    if dataset.type != "simplified":
        return {"lmax": -1.0, "muhat": None, "sigma_mu": None}
    args={"marginalize":marginalize}
    nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
    bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
    bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
    cov = dataset.globalInfo.covariance
    if type(nsig)==tuple:
        nsig = np.array(nsig)

    thirdMoment = None # Need to be implemented
    xsec = sum(nsig)/dataset.getLumi()

    statModel = get_multi_region_statistical_model(analysis=dataset.globalInfo.id,
                                                    signal=nsig,
                                                    observed=nobs,
                                                    covariance=cov,
                                                    nb=bg,
                                                    third_moment=thirdMoment,
                                                    delta_sys=deltas_rel,
                                                    xsection=xsec.asNumber(pb)
                                                    )
    expectedDict = {False:ExpectationType.observed,
                    True:ExpectationType.apriori,
                    "posteriori":ExpectationType.aposteriori}
    if expected not in expectedDict.keys():
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    config = statModel.backend.model.config()
    bounds = [(suggested[0]-200,suggested[1]+200) for suggested in config.suggested_bounds]
    if allowNegativeSignals:
        bounds[config.poi_index] = (config.minimum_poi, 100)
    else:
        bounds[config.poi_index] = (0, 100)
    muhat, lmax = statModel.maximize_likelihood(allow_negative_signal=allowNegativeSignals, expected=expectedDict[expected], return_nll=nll, par_bounds=bounds)
    lbsm = statModel.likelihood ( poi_test = 1., expected=expectedDict[expected], return_nll = nll, **args )
    lsm = statModel.likelihood ( poi_test = 0., expected=expectedDict[expected], return_nll = nll, **args )
    if lsm > lmax:
        logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lsm
        muhat = 0.0
    if lbsm > lmax:
        logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
        lmax = lbsm
        muhat = 1.0
    test_statistics = "q" if allowNegativeSignals else "qmutilde"
    sigma_mu = statModel.sigma_mu(poi_test=muhat,expected=expectedDict[expected],test_statistics=test_statistics)
    return {"muhat": muhat, "sigma_mu": sigma_mu, "lmax": lmax, "lbsm": lbsm, "lsm": lsm }


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

def getWSInfo(jsons):
    """
    Getting informations from the json files

    :param jsons: list of json instances.
    :return: wsInfo list of dictionaries (one dictionary for each json file) containing useful information about the json files.
        - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
        - :key otherRegions: list of strings indicating the path to the control and validation region channels
    """
    # Identifying the path to the SR and VR channels in the main workspace files
    wsInfo = []  # workspace specifications
    if not isinstance(jsons, list):
        logger.error("The `jsons` parameter must be of type list")
        return
    for ws in jsons:
        wsChannelsInfo = {}
        wsChannelsInfo["signalRegions"] = []
        wsChannelsInfo["otherRegions"] = []
        if not "channels" in ws.keys():
            logger.error(
                "Json file number {} is corrupted (channels are missing)".format(
                    jsons.index(ws)
                )
            )
            wsInfo = None
            return
        for i_ch, ch in enumerate(ws["channels"]):
            if ch["name"][:2] == "SR":  # if channel name starts with 'SR'
                wsChannelsInfo["signalRegions"].append(
                    {
                        "path": "/channels/"
                        + str(i_ch)
                        + "/samples/0",  # Path of the new sample to add (signal prediction)
                        "size": len(ch["samples"][0]["data"]),
                    }
                )  # Number of bins
            else:
                wsChannelsInfo["otherRegions"].append("/channels/" + str(i_ch))
        wsChannelsInfo["otherRegions"].sort(
            key=lambda path: path.split("/")[-1], reverse=True
        )  # Need to sort correctly the paths to the channels to be removed
        wsInfo.append(wsChannelsInfo)
    return wsInfo

def patchMaker(jsons, wsInfo, nsignals, includeCRs):
    """
    Method that creates the list of patches to be applied to the `jsons` workspaces, one for each region given the `nsignals` and the informations available in `wsInfo` and the content of the `jsons`

    :param jsons: list of json instances.
    :param wsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
    :param nsignals: list of list of signal yields (one list for each json file).
    :param includeCRs: if True leaves the pacth unchanged,
                       if False adds to the patch an operation that removes the CRs from the json files.

    NB: It seems we need to include the change of the "modifiers" in the patches as well

    :return: the list of patches, one for each workspace
    """
    if wsInfo == None:
        return None
    # Constructing the patches to be applied on the main workspace files
    patches = []
    for ws, info, subSig in zip(jsons, wsInfo, nsignals):
        patch = []
        for srInfo in info["signalRegions"]:
            nBins = srInfo["size"]
            operator = {}
            operator["op"] = "add"
            operator["path"] = srInfo["path"]
            value = {}
            value["data"] = subSig[:nBins]
            subSig = subSig[nBins:]
            value["modifiers"] = []
            value["modifiers"].append({"data": None, "type": "normfactor", "name": "mu_SIG"})
            value["modifiers"].append({"data": None, "type": "lumi", "name": "lumi"})
            value["name"] = "bsm"
            operator["value"] = value
            patch.append(operator)
        if includeCRs:
            logger.debug("keeping the CRs")
        else:
            for path in info["otherRegions"]:
                patch.append({"op": "remove", "path": path})
        patches.append(patch)
    return patches

def _getPatches(dataset, nsig):
    """
    :param nsig: list of signal yields (not relative).
    :param normalize: if true, normalize nsig
    :returns: the list of patches, one for each wkspace and the list of list of signals (one list for each json file).
    """
    # Getting the path to the json files
    datasets = [ds.getID() for ds in dataset.origdatasets]
    jsonFiles = [js for js in dataset.globalInfo.jsonFiles] #List of json files names (list of str)
    jsons = dataset.globalInfo.jsons.copy() #List of json files (list of dict)
    if not isinstance(jsons, list):
        logger.error("The `jsons` parameter must be of type list.")
        return None
    # Filtering the json files by looking at the available datasets
    listOfSRInJson = []
    for jsName in dataset.globalInfo.jsonFiles:
        if all([ds not in dataset.globalInfo.jsonFiles[jsName] for ds in datasets]):
            # No datasets found for this json combination
            jsIndex = jsonFiles.index(jsName)
            jsonFiles.pop(jsIndex)
            jsons.pop(jsIndex)
            continue
        if not all([ds in datasets for ds in dataset.globalInfo.jsonFiles[jsName]]):
            # Some SRs are missing for this json combination
            logger.error( "Wrong json definition in globalInfo.jsonFiles for json: %s" % jsName)
        listOfSRInJson += dataset.globalInfo.jsonFiles[jsName]
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

    if hasattr(dataset.globalInfo, "includeCRs"):
        includeCRs = dataset.globalInfo.includeCRs
    else:
        includeCRs = False

    wsInfo = getWSInfo(jsons)
    return patchMaker(jsons, wsInfo, nsignals, includeCRs), nsignals
